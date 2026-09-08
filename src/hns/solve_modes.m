function [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff] = solve_modes(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Solve the harmonic Navier-Stokes (HNS) system for all active modes.
%
%   Opt.lu_mode selects the linear-solver back-end (LU = LU factorization):
%     'lapack_band'    banded LU (compiled MEX routine), factor cached per mode (compressible)
%     'block_thomas'   dual-Schur block-Thomas (Paredes 2024), compressible only
%     'full'           cache every LU factor (fastest reuse, most memory)
%     'auto'           least-recently-used (LRU) cache bounded by Opt.lu_max_cache
%     'seq'            one factor at a time, discard after solve (low memory)
%     'none'           backslash each solve (no caching)
%   Opt.lu_cache='off' forces the cached back-ends to one-shot factor+solve.
%   Opt.parfor='on' distributes the mode solves (auto-set in init).
%   Opt.bt_apply selects how block_thomas applies the diagonal station block:
%     'lu' (default) stores its pivoted LU factors, 'inverse' stores G^-1.

%% Inputs
% Name     size              unit   explanation
% BF       (1,1) struct      [-]    base flow on its own grid
% Bound    (1,1) struct      [-]    boundary conditions and outflow-buffer coefficients
% StabGrid (1,1) struct      [-]    stability grid, coordinate metrics and base-flow fields on it
% Init     (1,1) struct      [-]    inflow initialization data
% Mode     (1,1) struct      [-]    mode bookkeeping: nf, active set RunJ and the interaction tables
% Opt      (1,1) struct      [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% Stab     (1,1) struct      [-]    disturbance specification: modal truncation, frequencies, inflow condition and amplitudes
% StabRes  (1,1) struct      [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% MLA1     nvar*ny*nx sparse [-]    mode-invariant operator block, zeroth order in omega and beta
% MLA2     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in omega
% MLA3     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in beta
% MLA4     nvar*ny*nx sparse [-]    mode-invariant operator block, second order in beta
% MLA5     nvar*ny*nx sparse [-]    mode-invariant operator block, mixed omega-beta term
% MFD_diff nvar*ny*nx sparse [-]    difference operator for the mean-flow-distortion block

%% Outputs
% Name     size              unit   explanation
% StabRes  (1,1) struct      [-]    solution: shape functions, modal amplitudes, wavenumbers and iteration history
% Opt      (1,1) struct      [-]    solver options; see the option index in docs/DeHNSSo_guide.md
% MLA1     nvar*ny*nx sparse [-]    mode-invariant operator block, zeroth order in omega and beta
% MLA2     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in omega
% MLA3     nvar*ny*nx sparse [-]    mode-invariant operator block, first order in beta
% MLA4     nvar*ny*nx sparse [-]    mode-invariant operator block, second order in beta
% MLA5     nvar*ny*nx sparse [-]    mode-invariant operator block, mixed omega-beta term
% MFD_diff nvar*ny*nx sparse [-]    difference operator for the mean-flow-distortion block

if ~isfield(Opt, 'lu_mode'); Opt.lu_mode = 'auto'; end
if ~isfield(Opt, 'parfor');  Opt.parfor  = 'off';  end
% Fill-reducing ordering for the sparse-LU back-ends: 'auto' (default) picks
% symrcm for the compressible operator and COLAMD for the incompressible one,
% which is the better choice on each (see do_lu). 'symrcm' or 'colamd' force one.
if ~isfield(Opt, 'lu_ordering') || isempty(Opt.lu_ordering); Opt.lu_ordering = 'auto'; end
% Opt.lu_scaling='on' keeps UMFPACK's row scaling by taking the five-output lu.
% The reduced-output forms below drop it, which costs accuracy on this operator
% because its rows differ in scale by many decades.
use_scaling = isfield(Opt,'lu_scaling') && strcmpi(Opt.lu_scaling,'on');
% Opt.bt_apply chooses how block_thomas applies the diagonal station block:
% 'inverse' (default) stores G^-1 as in Paredes et al. (2024), 'lu' keeps the
% pivoted factors instead. See block_thomas for the measured difference.
if ~isfield(Opt,'bt_apply') || isempty(Opt.bt_apply); Opt.bt_apply = 'lu'; end
bt_amode = Opt.bt_apply;

phys_modes = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2));

%% Phase 1: Build LHS matrices (once per run)
needs_beta = any(StabRes.betavec ~= 0);
if Opt.mmry == 1 && needs_beta && (isempty(MLA3) || isempty(MLA4))
    fprintf(' MLA3/MLA4 missing for 3D run, rebuilding LHS.\n');
    Opt.mmry = 0;
end
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible, 'on');
if Opt.mmry ~= 1
    Opt.mmry = 1;
    Opt.hr_in_MLA1 = false;
    Opt.hr_rows    = [];
    fprintf(' Constructing LHS matrices ...\n');
    tic
    if is_compr
        [MLAMFD, MLA1, MLA2, MLA3, MLA4, MLA5] = LHS_compressible(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
    else
        [MLAMFD, MLA1, MLA2, MLA3, MLA4]       = LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
        MLA5 = sparse(size(MLA1,1), size(MLA1,2));
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

%% Phase 1b: Resonator wall rows (also applied to a restored LHS)
% Helmholtz-resonator impedance wall boundary condition (BC), plus an optional
% free-slip u' over the mouth. Applied to MLA1 only, so the (0,0) mean-flow
% distortion (MFD) mode keeps no-slip while the fluctuations see the resonator
% (MFD_diff = MLAMFD - MLA1 restores no-slip for the MFD mode).
want_hr = isfield(Bound,'hr') && ~isempty(Bound.hr);
if ~Opt.hr_in_MLA1 && (want_hr || ~isempty(Opt.hr_rows))
    MLAMFD = MLA1 + MFD_diff;                        % invariant under the wall rows
    if ~isempty(Opt.hr_rows)
        MLA1(Opt.hr_rows, :) = MLAMFD(Opt.hr_rows, :);
        Opt.hr_rows = [];
    end
    if want_hr
        [MLA1, Opt.hr_rows] = apply_resonator_wall(MLA1, Bound.hr, BF, StabGrid, is_compr);
    end
    MFD_diff = MLAMFD - MLA1;
    Opt.hr_in_MLA1 = want_hr;   % Picard calls solve_modes per iteration; write once
    clear MLAMFD
end

%% Phase 2 & 3: Factorize and solve
nphys  = numel(phys_modes);
q_out  = zeros(nphys, size(StabRes.q, 2));
Stab_f = Stab.f;
% Zero the nonlinear-term (NLT) forcing on all Dirichlet BC rows (single-entry).
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
    % Sequential mode: factorize one matrix, solve all modes that use it,
    % discard.  Only one LU in memory at a time: minimal memory footprint.
    unique_keys = {};
    mode_groups = {};
    for k = 1:nphys
        j = phys_modes(k);
        key = ml_key(StabRes.omegavec(j), StabRes.betavec(j));
        idx = find(strcmp(unique_keys, key), 1);
        if isempty(idx)
            unique_keys{end+1} = key;
            mode_groups{end+1} = k;
        else
            mode_groups{idx}(end+1) = k;
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
        % Same factorization and ordering choice as the cached back-ends, so that
        % Opt.lu_ordering reaches this path too; lu_solve handles either layout.
        lu_s = do_lu(ML, pick_ordering(Opt, is_compr), use_scaling);

        % Iterative refinement (Opt.iterref passes, default 0 = off). The operator
        % is severely ill-conditioned, so plain LU leaves deeply sub-dominant modes
        % at ~kappa*eps; refinement recovers accuracy towards eps. Reuses the LU.
        niter = 0;
        if isfield(Opt,'iterref') && ~isempty(Opt.iterref); niter = Opt.iterref; end
        if niter > 0; MLr = ML; end          % refinement needs the unpermuted operator
        clear ML

        for ki = 1:numel(k_list)
            k = k_list(ki);
            j = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_f(j,:);
            f_j(R_j ~= 0) = 0;
            rhs = R_j + f_j.';
            x   = lu_solve(lu_s, rhs);
            for it = 1:niter                  % refine on the original ordering
                x = x + lu_solve(lu_s, rhs - MLr*x);
            end
            q_out(k,:) = x.';
        end
        clear lu_s MLr
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
            unique_keys{end+1} = key;
            mode_groups{end+1} = k;
        else
            mode_groups{idx}(end+1) = k;
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
            unique_keys{end+1} = key; mode_groups{end+1} = k;
        else
            mode_groups{idx}(end+1) = k;
        end
    end
    ntot = size(MLA1,1); b = round(ntot / StabGrid.nx);

    % Parallel pre-factorization (Opt.parfor='on', in-core caching). The block
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
            amode_p = bt_amode;                 % broadcast into the parfor body
            parfor ii = 1:numel(fi)
                MLp = assemble_ML(omgf(ii), betf(ii), MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
                btF_new{ii} = block_thomas('factor', MLp, b, np_, '', [], amode_p);
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
            Xsol = block_thomas('solve', ML, b, RHSm, Opt.bt_np, oocdir, [], bt_amode);
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
                Opt.LU_cache.(key) = struct('btF', block_thomas('factor', ML, b, Opt.bt_np, oocdir, [], bt_amode));
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

    % Keys needed by the solve loop this iteration: 'auto' eviction must not
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
                    Opt.LU_cache.(ML_key) = do_lu(ML, pick_ordering(Opt, is_compr), use_scaling);
                    fprintf(' (%.1f s, L+U: %.0f MB)', toc(tic_lu), ...
                        (nnz(Opt.LU_cache.(ML_key).L)+nnz(Opt.LU_cache.(ML_key).U))*16/1e6);

                case 'auto'
                    if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = max(8, nphys); end
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
                            % Cache saturated with needed keys: must grow to avoid
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
                    Opt.LU_cache.(ML_key) = do_lu(ML, pick_ordering(Opt, is_compr), use_scaling);
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
        Stab_fk = Stab_f(phys_modes, :);   % row per worker: sliced, not broadcast
        parfor k = 1:nphys
            j   = phys_modes(k);
            R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
            f_j = Stab_fk(k,:);
            f_j(R_j ~= 0) = 0;
            q_out(k,:) = lu_solve(LU_cell{k}, R_j + f_j.').';
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
if isstruct(lu_s) && isfield(lu_s, 'D')
    % Five-output lu: P*(D\ML)*Q = L*U  =>  x = Q*(U\(L\(P*(D\rhs))))
    if isfield(lu_s, 'r')
        y = lu_s.Q * (lu_s.U \ (lu_s.L \ (lu_s.P * (lu_s.D \ rhs(lu_s.r)))));
        x = zeros(size(rhs));
        x(lu_s.r) = y;
    else
        x = lu_s.Q * (lu_s.U \ (lu_s.L \ (lu_s.P * (lu_s.D \ rhs))));
    end
elseif isstruct(lu_s) && isfield(lu_s, 'q')
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

function s = do_lu(ML, use_symrcm, use_scaling)
%DO_LU  LU with a fill-reducing ordering.
%   use_symrcm=true: symrcm + 3-output lu. The operator is a chain of station
%     blocks, so it is banded and bandwidth reduction is what limits the fill;
%     needed on the COMPRESSIBLE operator, where COLAMD fills catastrophically.
%   use_symrcm=false: 4-output lu (COLAMD with column pivoting), which gives both
%     less fill and a faster factorization on the INCOMPRESSIBLE operator.
%   Both are handled by lu_solve. The caller decides; solve_modes picks by
%   operator unless Opt.lu_ordering overrides it.
%
%   Measured fill, nnz(L+U)/nnz(M), at 40/60/80 stations:
%     incompressible  COLAMD 6.7/6.8/6.8   symrcm 7.1/7.3/7.3
%     compressible    COLAMD 25/43/50      symrcm 36.4/37.1/37.7
%   symrcm stays flat in Nx while COLAMD grows, crossing near Nx=55; at Nx=200
%   COLAMD gave nearly double the fill. AMD and nested dissection are far worse
%   on both (up to 79x and 174x), so neither is offered.
if nargin < 2 || isempty(use_symrcm);  use_symrcm  = false; end
if nargin < 3 || isempty(use_scaling); use_scaling = false; end
if use_scaling
    % Five-output lu: P*(D\ML)*Q = L*U, so the row scaling D is kept.
    if use_symrcm
        r = symrcm(ML);
        [L, U, P, Q, D] = lu(ML(r,r));
        s = struct('L', L, 'U', U, 'P', P, 'Q', Q, 'D', D, 'r', r);
    else
        [L, U, P, Q, D] = lu(ML);
        s = struct('L', L, 'U', U, 'P', P, 'Q', Q, 'D', D);
    end
elseif use_symrcm
    r = symrcm(ML);
    [L, U, p] = lu(ML(r,r), 'vector');
    s = struct('L', L, 'U', U, 'p', p, 'r', r);
else
    [L, U, p, q] = lu(ML, 'vector');
    s = struct('L', L, 'U', U, 'p', p, 'q', q);
end
end

function use_symrcm = pick_ordering(Opt, is_compr)
%PICK_ORDERING  Resolve Opt.lu_ordering to the flag do_lu expects.
%   'auto' (default) follows the operator: symrcm for the compressible one,
%   COLAMD for the incompressible one. 'symrcm' and 'colamd' force the choice,
%   which is worth doing on small compressible cases, where COLAMD is both
%   leaner and faster below roughly Nx = 55.
switch lower(Opt.lu_ordering)
    case 'auto';   use_symrcm = is_compr;
    case 'symrcm'; use_symrcm = true;
    case 'colamd'; use_symrcm = false;
    otherwise
        error('solve_modes:lu_ordering', ...
            'Opt.lu_ordering must be ''auto'', ''symrcm'' or ''colamd'' (got ''%s'')', ...
            Opt.lu_ordering);
end
end
