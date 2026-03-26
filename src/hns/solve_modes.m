function [StabRes, Opt, MLA1, MLA2, MFD_diff] = solve_modes(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MFD_diff)
%SOLVE_MODES  Solve the HNS system for all active modes.
%
%   LU caching controlled by Opt.lu_mode:
%     'full' — cache all LU factorisations (fastest, most memory)
%     'auto' — cache up to Opt.lu_max_cache, evict oldest (balanced)
%     'none' — use backslash each time, no caching (slowest, least memory)

if ~isfield(Opt, 'lu_mode'); Opt.lu_mode = 'full'; end

for j = [Mode.RunJ(Mode.RunJ >= round(Mode.nf/2))]
    omega = StabRes.omegavec(j);
    beta  = StabRes.betavec(j);

    % Build LHS matrices (once)
    if Opt.mmry ~= 1
        Opt.mmry = 1;
        fprintf(' Constructing LHS matrices ...\n');
        tic
        [MLAMFD, MLA1, MLA2] = LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
        MFD_diff = MLAMFD - MLA1;
        clear MLAMFD
        fprintf(' LHS done in %.1f s\n', toc);
    end

    fprintf(' Solving (%d,%d)', Stab.Mvec(j), Stab.Nvec(j));

    % Cache key
    if omega == 0 && beta == 0
        ML_key = 'MFD';
    else
        ML_key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
    end

    if ~isfield(Opt, 'LU_cache'); Opt.LU_cache = struct(); end
    need_factor = ~isfield(Opt.LU_cache, ML_key);

    if need_factor
        % Assemble ML
        if omega == 0 && beta == 0
            ML = MLA1 + MFD_diff;
        elseif omega == 0
            ML = MLA1 + MLA3*beta + MLA4*beta^2;
        elseif beta == 0
            ML = MLA1 + MLA2*omega;
        else
            ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2;
        end

        switch Opt.lu_mode
            case 'full'
                fprintf(' [LU: cache all]');
                Opt.LU_cache.(ML_key) = decomposition(ML, 'lu');

            case 'auto'
                if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = 3; end
                if ~isfield(Opt, 'lu_order'); Opt.lu_order = {}; Opt.lu_order_label = {}; end
                n_cached = numel(fieldnames(Opt.LU_cache));
                if n_cached >= Opt.lu_max_cache
                    old_key = Opt.lu_order{1};
                    old_mode = Opt.lu_order_label{1};
                    Opt.LU_cache = rmfield(Opt.LU_cache, old_key);
                    Opt.lu_order(1) = [];
                    Opt.lu_order_label(1) = [];
                    fprintf(' [LU: replacing %s, %d cached]', old_mode, Opt.lu_max_cache);
                else
                    fprintf(' [LU: %d/%d cached]', n_cached+1, Opt.lu_max_cache);
                end
                Opt.LU_cache.(ML_key) = decomposition(ML, 'lu');
                Opt.lu_order{end+1} = ML_key;
                Opt.lu_order_label{end+1} = sprintf('(%d,%d)', Stab.Mvec(j), Stab.Nvec(j));

            case 'none'
                fprintf(' [direct]');
                Opt.LU_cache.(ML_key) = ML;
        end
        clear ML
    end

    % Construct RHS
    [StabRes.R] = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
    Stab.f(j, StabRes.R ~= 0) = 0;

    % Solve
    rhs = StabRes.R + Stab.f(j,:).';
    StabRes.q(j,:) = (Opt.LU_cache.(ML_key) \ rhs).';
    StabRes.q(Mode.nf+1-j,:) = conj(StabRes.q(j,:));

    fprintf(' ... %.1fs\n', toc)
end

end
