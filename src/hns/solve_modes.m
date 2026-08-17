function [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff] = solve_modes(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
%SOLVE_MODES  Solve the HNS system for all active modes.
%
%   Opt.lu_mode selects the linear-solver back-end:
%     'lapack_band'  — banded LU MEX, factor cached per mode (compressible)
%     'block_thomas' — dual-Schur block-Thomas (Paredes 2024), compressible only
%     'full'         — cache every LU factor (fastest reuse, most memory)
%     'auto'         — LRU cache bounded by Opt.lu_max_cache
%     'seq'          — one factor at a time, discard after solve (low memory)
%     'none'         — backslash each solve (no caching)
%   Opt.lu_cache='off' forces the cached back-ends to one-shot factor+solve.
%   Opt.parfor='on' distributes the mode solves (auto-set in init).

if ~isfield(Opt, 'lu_mode'); Opt.lu_mode = 'full'; end
if ~isfield(Opt, 'parfor');  Opt.parfor  = 'off';  end

phys_modes = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2));

%% Phase 1: Build LHS matrices (once per run)
needs_beta = any(StabRes.betavec ~= 0);
if Opt.mmry == 1 && needs_beta && (isempty(MLA3) || isempty(MLA4))
    fprintf(' MLA3/MLA4 missing for 3D run — rebuilding LHS.\n');
    Opt.mmry = 0;
end
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible, 'on');
if Opt.mmry ~= 1
    Opt.mmry = 1;
    fprintf(' Constructing LHS matrices ...\n');
    tic
    if is_compr
        [MLAMFD, MLA1, MLA2, MLA3, MLA4, MLA5] = LHS_compressible(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
    else
        [MLAMFD, MLA1, MLA2, MLA3, MLA4]       = LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
        MLA5 = sparse(size(MLA1,1), size(MLA1,2));
    end
    % Helmholtz-resonator impedance wall BC (+ optional free-slip u' over the mouth).
    % Applied to MLA1 only, so the (0,0) MFD keeps no-slip while the fluctuations
    % see the resonator (MFD_diff = MLAMFD - MLA1 restores no-slip for the MFD mode).
    if isfield(Opt,'hr') && ~isempty(Opt.hr)
        MLA1 = apply_resonator_wall(MLA1, Opt.hr, BF, StabGrid, is_compr);
    end
    MFD_diff = MLAMFD - MLA1;
    clear MLAMFD
    fprintf(' LHS done in %.1f s\n', toc);
    if ~isempty(getenv('MFD_OPDUMP'))
        ML_mfd = MLA1 + MFD_diff;                       % MLAMFD (omega=beta=0)
        mfd_j  = find(Stab.Mvec==0 & Stab.Nvec==0, 1);
        R_mfd  = RHS(Bound, StabGrid, Mode, Opt, StabRes, mfd_j);
        f_mfd  = Stab.f(mfd_j,:);
        q_mfd  = StabRes.q(mfd_j,:);
        Xg = StabGrid.x(:).'; Yg = StabGrid.y(:); nyg = numel(Yg); nxg = numel(Xg);
        save('/tmp/mfd_op.mat','ML_mfd','R_mfd','f_mfd','q_mfd','Xg','Yg','nyg','nxg','mfd_j','-v7.3');
        fprintf(' [MFD_OPDUMP] saved MFD operator+forcing -> /tmp/mfd_op.mat\n');
    end
end

%% Phase 2 & 3: Factorise and solve
nphys  = numel(phys_modes);
q_out  = zeros(nphys, size(StabRes.q, 2));
Stab_f = Stab.f;
% Zero the NLT forcing on all Dirichlet BC rows (single-entry). The
% f_j(R_j~=0)=0 masks catch only inhomogeneous BCs, so homogeneous wall rows
% would otherwise inherit forcing (e.g. spurious (0,0) wall transpiration).
bc_dir     = full(sum(MLA1 ~= 0, 2)) == 1;
bc_dir_mfd = full(sum((MLA1 + MFD_diff) ~= 0, 2)) == 1;
for jj = 1:size(Stab_f, 1)
    if Stab.Mvec(jj) == 0 && Stab.Nvec(jj) == 0
        Stab_f(jj, bc_dir_mfd) = 0;
    else
        Stab_f(jj, bc_dir) = 0;
    end
end
Mvec   = Stab.Mvec;
Nvec   = Stab.Nvec;

if strcmpi(Opt.lu_mode, 'seq')
    % Sequential mode: factorise one matrix, solve all modes that use it,
    % discard.  Only one LU in memory at a time — minimal memory footprint.
    unique_keys = {};
    mode_groups = {};
    for k = 1:nphys
        j = phys_modes(k);
        key = ml_key(StabRes.omegavec(j), StabRes.betavec(j));
        idx = find(strcmp(unique_keys, key), 1);
        if isempty(idx)
            unique_keys{end+1} = key; %#ok<AGROW>
            mode_groups{end+1} = k;   %#ok<AGROW>
        else
            mode_groups{idx}(end+1) = k; %#ok<AGROW>
        end
    end

    for g = 1:numel(unique_keys)
        k_list = mode_groups{g};
        j0 = phys_modes(k_list(1));
        omega = StabRes.omegavec(j0);
        beta  = StabRes.betavec(j0);
        fprintf(' LU+solve (%d,%d) [seq, %d mode(s)] ...', Mvec(j0), Nvec(j0), numel(k_list));
        tic_g = tic;

        ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
        r  = symrcm(ML);
        [L, U, p] = lu(ML(r,r), 'vector');

        % Iterative refinement (Opt.iterref passes, default 0 = off). The operator
        % is severely ill-conditioned, so plain LU leaves deeply sub-dominant modes
        % at ~kappa*eps; refinement recovers accuracy towards eps. Reuses the LU.
        niter = 0;
        if isfield(Opt,'iterref') && ~isempty(Opt.iterref); niter = Opt.iterref; end
        if niter > 0; Arr = ML(r,r); end
        clear ML

        for ki = 1:numel(k_list)
            k = k_list(ki);
            j = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            rhs   = R_j + f_j.';
            rhs_r = rhs(r);
            x_r   = U \ (L \ rhs_r(p));
            for it = 1:niter
                res = rhs_r - Arr*x_r;
                x_r = x_r + U \ (L \ res(p));
            end
            x     = zeros(size(rhs));
            x(r)  = x_r;
            q_out(k,:) = x.';
        end
        clear L U p r
        fprintf(' done (%.1f s)\n', toc(tic_g));
    end

elseif strcmpi(Opt.lu_mode, 'lapack_band')
    % Banded LU with per-mode factor caching: factor once (dgbsv_factor), reuse
    % across Picard iters (dgbsv_solve_cached). Cache capped by Opt.lu_max_cache
    % (LRU, skips keys needed this iter); persists across runs via Opt.rerun.
    if exist(['dgbsv_factor.' mexext],'file') ~= 3 || ...
       exist(['dgbsv_solve_cached.' mexext],'file') ~= 3
        error('solve_modes:no_mex', ...
              'banded-LU MEX not built. Run src/linsolve/lapack_band/mexbuild_dgbsv.m first.');
    end
    if ~isfield(Opt, 'LU_cache'); Opt.LU_cache = struct(); end
    if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = max(8, nphys); end  % default: hold all modes
    % Banded factor caching on/off (matches the 'on'/'off' Opt convention).
    % 'off' (or lu_max_cache<=0) => one-shot factor+solve, one band resident.
    band_nocache = Opt.lu_max_cache <= 0 || ...
                   (isfield(Opt,'lu_cache') && strcmpi(Opt.lu_cache,'off'));

    unique_keys = {};
    mode_groups = {};
    for k = 1:nphys
        j = phys_modes(k);
        key = ml_key(StabRes.omegavec(j), StabRes.betavec(j));
        idx = find(strcmp(unique_keys, key), 1);
        if isempty(idx)
            unique_keys{end+1} = key; %#ok<AGROW>
            mode_groups{end+1} = k;   %#ok<AGROW>
        else
            mode_groups{idx}(end+1) = k; %#ok<AGROW>
        end
    end

    for g = 1:numel(unique_keys)
        key   = unique_keys{g};
        k_list = mode_groups{g};
        j0 = phys_modes(k_list(1));
        omega = StabRes.omegavec(j0);
        beta  = StabRes.betavec(j0);

        if band_nocache
            % No-cache (low-memory): one-shot factor+solve, one band at a time
            % (~1 band resident instead of nphys). Reproduces the pre-caching
            % footprint for memory-constrained machines.
            fprintf(' Band-LU solve (%d,%d) [lapack_band, no-cache]\n', Mvec(j0), Nvec(j0));
            ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
            niter = 0;
            if isfield(Opt,'iterref') && ~isempty(Opt.iterref); niter = Opt.iterref; end
            if niter > 0
                % Iterative refinement: factor once, refine each RHS with cheap
                % cached re-solves (one band resident, so memory is unchanged).
                [ABf, IPIV, KLb, KUb] = dgbsv_factor(ML);
                for ki = 1:numel(k_list)
                    k = k_list(ki); j = phys_modes(k);
                    R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
                    f_j = Stab_f(j,:); f_j(R_j ~= 0) = 0;
                    rhs = full(R_j + f_j.');
                    x = dgbsv_solve_cached(ABf, IPIV, KLb, KUb, rhs);
                    for it = 1:niter
                        x = x + dgbsv_solve_cached(ABf, IPIV, KLb, KUb, rhs - ML*x);
                    end
                    q_out(k,:) = x.';
                end
                clear ABf IPIV
            else
                for ki = 1:numel(k_list)
                    k = k_list(ki); j = phys_modes(k);
                    R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
                    f_j = Stab_f(j,:); f_j(R_j ~= 0) = 0;
                    rhs = full(R_j + f_j.');
                    q_out(k,:) = dgbsv_solve(ML, rhs).';
                end
            end
            clear ML
            continue
        end

        cached = isfield(Opt.LU_cache, key) && isstruct(Opt.LU_cache.(key)) ...
                 && isfield(Opt.LU_cache.(key), 'ABf');
        if ~cached
            % Evict (LRU) if at capacity, skipping keys needed this iteration.
            ck = fieldnames(Opt.LU_cache);
            if numel(ck) >= Opt.lu_max_cache
                evictable = ck(~ismember(ck, unique_keys));
                if ~isempty(evictable)
                    Opt.LU_cache = rmfield(Opt.LU_cache, evictable{1});
                end
            end
            fprintf(' Band-LU factor (%d,%d) [lapack_band] ...', Mvec(j0), Nvec(j0));
            tic_g = tic;
            ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
            [ABf, IPIV, KLb, KUb] = dgbsv_factor(ML);
            clear ML
            Opt.LU_cache.(key) = struct('ABf', ABf, 'IPIV', IPIV, 'KL', KLb, 'KU', KUb);
            fprintf(' done (%.1f s)\n', toc(tic_g));
        else
            fprintf(' Band-LU solve (%d,%d) [cached]\n', Mvec(j0), Nvec(j0));
        end

        C = Opt.LU_cache.(key);
        niter = 0;
        if isfield(Opt,'iterref') && ~isempty(Opt.iterref); niter = Opt.iterref; end
        if niter > 0   % re-assemble ML (sparse, cheap) for the refinement residual
            MLr = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
        end
        for ki = 1:numel(k_list)
            k = k_list(ki);
            j = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            rhs = full(R_j + f_j.');
            x = dgbsv_solve_cached(C.ABf, C.IPIV, C.KL, C.KU, rhs);
            for it = 1:niter
                x = x + dgbsv_solve_cached(C.ABf, C.IPIV, C.KL, C.KU, rhs - MLr*x);
            end
            q_out(k,:) = x.';
        end
        if niter > 0; clear MLr; end
    end

elseif strcmpi(Opt.lu_mode, 'block_thomas')
    % Dual-Schur block-Thomas (Paredes 2024), partitioned by Opt.bt_np (default 2).
    % Like the banded path, factor once per (omega,beta) and reuse via
    % block_thomas('factor'/'apply'); cached in Opt.LU_cache (LRU, Opt.rerun).
    %   Opt.lu_cache='off' => one-shot factor+solve.
    %   Opt.bt_ooc='on'    => out-of-core (factor streamed to disk, minimal RAM).
    if ~is_compr
        error('solve_modes:bt_incompressible', ...
            ['lu_mode=''block_thomas'' is for the COMPRESSIBLE operator only. The ' ...
             'incompressible saddle-point (zero-diagonal continuity rows) makes the ' ...
             'unpivoted block factorization near-singular and the result SILENTLY ' ...
             'WRONG (verified ~2x error). Use lu_mode=''lapack_band'', ''seq'', or ''auto''.']);
    end
    if ~isfield(Opt, 'bt_np') || isempty(Opt.bt_np); Opt.bt_np = 2; end
    if ~isfield(Opt, 'LU_cache'); Opt.LU_cache = struct(); end
    if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = max(8, nphys); end
    bt_ooc     = isfield(Opt,'bt_ooc') && strcmpi(Opt.bt_ooc,'on');
    bt_nocache = Opt.lu_max_cache <= 0 || ...
                 (isfield(Opt,'lu_cache') && strcmpi(Opt.lu_cache,'off'));
    unique_keys = {}; mode_groups = {};
    for k = 1:nphys
        j   = phys_modes(k);
        key = ml_key(StabRes.omegavec(j), StabRes.betavec(j));
        idx = find(strcmp(unique_keys, key), 1);
        if isempty(idx)
            unique_keys{end+1} = key; mode_groups{end+1} = k; %#ok<AGROW>
        else
            mode_groups{idx}(end+1) = k; %#ok<AGROW>
        end
    end
    ntot = size(MLA1,1); b = round(ntot / StabGrid.nx);

    % Parallel pre-factorisation (Opt.parfor='on', in-core caching). The block
    % factor is the expensive step (the apply is cheap), so factor every
    % not-yet-cached group CONCURRENTLY here; the serial loop below then finds
    % them cached and only applies. Needs a Threads parpool so MLA1..5 are shared
    % (not copied). Skipped for no-cache / OOC, which factor one-shot per solve.
    if strcmpi(Opt.parfor,'on') && ~bt_nocache && ~bt_ooc
        ngp = numel(unique_keys);
        need = false(ngp,1); omg = zeros(ngp,1); bet = zeros(ngp,1);
        for gp = 1:ngp
            j0g = phys_modes(mode_groups{gp}(1));
            omg(gp) = StabRes.omegavec(j0g); bet(gp) = StabRes.betavec(j0g);
            need(gp) = ~(isfield(Opt.LU_cache, unique_keys{gp}) ...
                         && isstruct(Opt.LU_cache.(unique_keys{gp})) ...
                         && isfield(Opt.LU_cache.(unique_keys{gp}), 'btF'));
        end
        fi = find(need);
        if ~isempty(fi)
            btF_new = cell(numel(fi),1); np_ = Opt.bt_np;
            omgf = omg(fi); betf = bet(fi);
            fprintf(' Block-Thomas: parallel factor of %d group(s) [np=%d] ...\n', numel(fi), np_);
            tf = tic;
            parfor ii = 1:numel(fi)
                MLp = assemble_ML(omgf(ii), betf(ii), MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
                btF_new{ii} = block_thomas('factor', MLp, b, np_, '');
            end
            fprintf(' Block-Thomas parallel factor done (%.1f s)\n', toc(tf));
            for ii = 1:numel(fi)
                Opt.LU_cache.(unique_keys{fi(ii)}) = struct('btF', btF_new{ii});
                btF_new{ii} = [];   % free as we move it into the cache (avoid a double-hold)
            end
        end
    end

    for g = 1:numel(unique_keys)
        key    = unique_keys{g};
        k_list = mode_groups{g};
        j0     = phys_modes(k_list(1));
        omega  = StabRes.omegavec(j0); beta = StabRes.betavec(j0);
        RHSm = zeros(ntot, numel(k_list));
        for ki = 1:numel(k_list)
            j   = phys_modes(k_list(ki));
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:); f_j(R_j ~= 0) = 0;
            RHSm(:,ki) = R_j + f_j.';
        end

        if bt_nocache
            oocstr = ''; if bt_ooc; oocstr = ', ooc'; end
            fprintf(' Block-Thomas solve (%d,%d) [np=%d, %d mode(s)%s] ...', ...
                    Mvec(j0), Nvec(j0), Opt.bt_np, numel(k_list), oocstr);
            tic_g = tic;
            ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
            oocdir = '';
            if bt_ooc; oocdir = fullfile(tempdir, sprintf('btooc_g%d', g)); end
            Xsol = block_thomas('solve', ML, b, RHSm, Opt.bt_np, oocdir);
            clear ML
            fprintf(' done (%.1f s)\n', toc(tic_g));
        else
            cached = isfield(Opt.LU_cache, key) && isstruct(Opt.LU_cache.(key)) ...
                     && isfield(Opt.LU_cache.(key), 'btF');
            if ~cached
                ck = fieldnames(Opt.LU_cache);
                if numel(ck) >= Opt.lu_max_cache
                    evictable = ck(~ismember(ck, unique_keys));
                    if ~isempty(evictable)
                        ev = evictable{1};
                        if isfield(Opt.LU_cache.(ev),'btF'); block_thomas('free', Opt.LU_cache.(ev).btF); end
                        Opt.LU_cache = rmfield(Opt.LU_cache, ev);
                    end
                end
                oocstr = ''; if bt_ooc; oocstr = ', ooc'; end
                fprintf(' Block-Thomas factor (%d,%d) [np=%d%s] ...', Mvec(j0), Nvec(j0), Opt.bt_np, oocstr);
                tic_g = tic;
                ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
                oocdir = '';
                if bt_ooc; oocdir = [tempname '_btcache']; block_thomas('free', oocdir,'reg'); end
                Opt.LU_cache.(key) = struct('btF', block_thomas('factor', ML, b, Opt.bt_np, oocdir));
                clear ML
                fprintf(' done (%.1f s)\n', toc(tic_g));
            else
                fprintf(' Block-Thomas solve (%d,%d) [cached]\n', Mvec(j0), Nvec(j0));
            end
            Xsol = block_thomas('apply', Opt.LU_cache.(key).btF, RHSm);
        end

        for ki = 1:numel(k_list)
            q_out(k_list(ki),:) = Xsol(:,ki).';
        end
        clear RHSm Xsol
    end

else
    % Cache-based modes: 'full', 'auto', 'none'
    if ~isfield(Opt, 'LU_cache'); Opt.LU_cache = struct(); end

    % Keys needed by the solve loop this iteration — 'auto' eviction must not
    % remove any of these (otherwise the solve below reads a missing entry).
    needed_keys = cell(1, nphys);
    for k = 1:nphys
        jk = phys_modes(k);
        needed_keys{k} = ml_key(StabRes.omegavec(jk), StabRes.betavec(jk));
    end
    needed_set = unique(needed_keys);

    for j = phys_modes
        omega   = StabRes.omegavec(j);
        beta    = StabRes.betavec(j);
        ML_key  = ml_key(omega, beta);

        if ~isfield(Opt.LU_cache, ML_key)
            ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);

            switch Opt.lu_mode
                case 'full'
                    fprintf(' LU (%d,%d) [cache all] ...', Stab.Mvec(j), Stab.Nvec(j));
                    tic_lu = tic;
                    Opt.LU_cache.(ML_key) = do_lu(ML, is_compr);
                    fprintf(' (%.1f s, L+U: %.0f MB)', toc(tic_lu), ...
                        (nnz(Opt.LU_cache.(ML_key).L)+nnz(Opt.LU_cache.(ML_key).U))*16/1e6);

                case 'auto'
                    if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = 3; end
                    if ~isfield(Opt, 'lu_order');     Opt.lu_order = {}; Opt.lu_order_label = {}; end
                    n_cached = numel(fieldnames(Opt.LU_cache));
                    if n_cached >= Opt.lu_max_cache
                        % Evict the oldest cached key that's NOT needed this iteration.
                        evict_idx = [];
                        for ii = 1:numel(Opt.lu_order)
                            if ~any(strcmp(Opt.lu_order{ii}, needed_set))
                                evict_idx = ii;
                                break;
                            end
                        end
                        if ~isempty(evict_idx)
                            old_key  = Opt.lu_order{evict_idx};
                            old_mode = Opt.lu_order_label{evict_idx};
                            Opt.LU_cache = rmfield(Opt.LU_cache, old_key);
                            Opt.lu_order(evict_idx) = [];
                            Opt.lu_order_label(evict_idx) = [];
                            fprintf(' LU (%d,%d) [replacing %s, %d cached] ...', Stab.Mvec(j), Stab.Nvec(j), old_mode, Opt.lu_max_cache);
                        else
                            % Cache saturated with needed keys — must grow to avoid
                            % thrash. lu_max_cache is a soft limit in this case.
                            if ~isfield(Opt,'lu_grew_warned') || ~Opt.lu_grew_warned
                                warning(['LU cache grown from lu_max_cache=%d to %d ' ...
                                         'to hold all %d active modes; set lu_max_cache ' ...
                                         'accordingly or use lu_mode=''seq''/''lapack_band'' ' ...
                                         'for a bounded footprint.'], ...
                                        Opt.lu_max_cache, n_cached+1, nphys);
                                Opt.lu_grew_warned = true;
                            end
                            fprintf(' LU (%d,%d) [%d cached, grown] ...', Stab.Mvec(j), Stab.Nvec(j), n_cached+1);
                        end
                    else
                        fprintf(' LU (%d,%d) [%d/%d cached] ...', Stab.Mvec(j), Stab.Nvec(j), n_cached+1, Opt.lu_max_cache);
                    end
                    Opt.LU_cache.(ML_key) = do_lu(ML, is_compr);   % ordering picked for the operator (see do_lu)
                    Opt.lu_order{end+1}       = ML_key;
                    Opt.lu_order_label{end+1} = sprintf('(%d,%d)', Stab.Mvec(j), Stab.Nvec(j));

                case 'none'
                    fprintf(' LU (%d,%d) [direct] ...', Stab.Mvec(j), Stab.Nvec(j));
                    Opt.LU_cache.(ML_key) = ML;
            end
            fprintf(' done\n');
            clear ML
        end
    end

    % Free MLA matrices once all unique modes are cached
    n_unique = numel(fieldnames(Opt.LU_cache));
    n_needed = numel(unique(arrayfun(@(j) ml_key(StabRes.omegavec(j), StabRes.betavec(j)), ...
        1:Mode.nf, 'UniformOutput', false)));
    if n_unique >= n_needed
        MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MLA5 = []; MFD_diff = [];
    end

    % Solve
    LU_cell = cell(nphys, 1);
    for k = 1:nphys
        j = phys_modes(k);
        LU_cell{k} = Opt.LU_cache.(ml_key(StabRes.omegavec(j), StabRes.betavec(j)));
    end

    if strcmpi(Opt.parfor, 'on')
        fprintf(' Solving %d modes in parallel ...\n', nphys);
        tic_solve = tic;
        parfor k = 1:nphys
            j   = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            q_out(k,:) = lu_solve(LU_cell{k}, R_j + f_j.').';   %#ok<PFBNS>
        end
        fprintf(' Solve done in %.1f s\n', toc(tic_solve));
    else
        for k = 1:nphys
            j = phys_modes(k);
            fprintf(' Solving (%d,%d)', Mvec(j), Nvec(j));
            tic_j = tic;
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            q_out(k,:) = lu_solve(LU_cell{k}, R_j + f_j.').';
            fprintf(' ... %.1fs\n', toc(tic_j));
        end
    end
end


%% Collect results (including conjugate modes)
for k = 1:nphys
    j = phys_modes(k);
    StabRes.q(j,:)              = q_out(k,:);
    StabRes.q(Mode.nf+1-j,:)   = conj(q_out(k,:));
end


if ~isempty(getenv('MFD_QDUMP'))
    mfd_j = find(Stab.Mvec==0 & Stab.Nvec==0, 1);
    kmfd  = find(phys_modes==mfd_j, 1);
    if ~isempty(kmfd)
        q_mfd_conv = q_out(kmfd,:);                     % CONVERGED MFD field (operator-ordered)
        Xg = StabGrid.x(:).'; Yg = StabGrid.y(:); nyg = numel(Yg); nxg = numel(Xg);
        save('/tmp/mfd_qconv.mat','q_mfd_conv','Xg','Yg','nyg','nxg','mfd_j','-v7.3');  % overwritten each Picard iter -> last = converged
    end
end

end

%% Local helpers

function key = ml_key(omega, beta)
if omega == 0 && beta == 0
    key = 'MFD';
else
    key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
end
end

function x = lu_solve(lu_s, rhs)
%LU_SOLVE  Solve using cached LU.
if isstruct(lu_s) && isfield(lu_s, 'q')
    % [L,U,p,q] = lu(ML,'vector'):  ML(p,q) = L*U  =>  x(q) = U\(L\rhs(p))
    x = zeros(size(rhs));
    x(lu_s.q) = lu_s.U \ (lu_s.L \ rhs(lu_s.p));
elseif isstruct(lu_s) && isfield(lu_s, 'dML')
    % symrcm + decomposition
    x = zeros(size(rhs));
    x(lu_s.r) = lu_s.dML \ rhs(lu_s.r);
elseif isstruct(lu_s) && isfield(lu_s, 'L')
    % symrcm + raw L,U,p
    rhs_r = rhs(lu_s.r);
    x = zeros(size(rhs));
    x(lu_s.r) = lu_s.U \ (lu_s.L \ rhs_r(lu_s.p));
else
    % decomposition object or raw sparse matrix
    x = lu_s \ rhs;
end
end

function s = do_lu(ML, use_symrcm)
%DO_LU  LU with the fill-reducing ordering chosen for the operator.
%   use_symrcm=true: symrcm + 3-output lu — needed on the COMPRESSIBLE operator
%     (COLAMD gives catastrophic fill there). false: 4-output lu (COLAMD), better
%     on the INCOMPRESSIBLE operator. Both handled by lu_solve.
if nargin < 2 || isempty(use_symrcm); use_symrcm = false; end
if use_symrcm
    r = symrcm(ML);
    [L, U, p] = lu(ML(r,r), 'vector');
    s = struct('L', L, 'U', U, 'p', p, 'r', r);
else
    [L, U, p, q] = lu(ML, 'vector');
    s = struct('L', L, 'U', U, 'p', p, 'q', q);
end
end

function ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
% MLA5 = iβ·∂_x cross-viscous block (compressible only; zero sparse otherwise).
if omega == 0 && beta == 0
    ML = MLA1 + MFD_diff;
elseif omega == 0
    ML = MLA1 + MLA3*beta + MLA4*beta^2 + MLA5*beta;
elseif beta == 0
    ML = MLA1 + MLA2*omega;
else
    ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2 + MLA5*beta;
end
end
