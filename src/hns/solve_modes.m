function [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MFD_diff] = solve_modes(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MFD_diff)
%SOLVE_MODES  Solve the HNS system for all active modes.
%
%   LU caching controlled by Opt.lu_mode:
%     'full' — cache all LU factorisations (fastest, most memory)
%     'auto' — cache up to Opt.lu_max_cache, evict oldest (balanced)
%     'none' — use backslash each time, no caching (slowest, least memory)
%
%   Parallel solve controlled by Opt.parfor (set automatically in init):
%     true  — parfor over modes (requires Parallel Computing Toolbox)
%     false — standard for loop
%
%   GPU solve controlled by Opt.gpu (set automatically in init):
%     true  — transfer ML to GPU, solve with gpuArray \ (requires PCT + NVIDIA GPU)
%             Note: GPU does not support decomposition caching; ML is stored as
%             gpuArray and refactorised each solve. Beneficial for large ny (>>80).
%     false — CPU solve (default)

if ~isfield(Opt, 'lu_mode'); Opt.lu_mode = 'full'; end
if ~isfield(Opt, 'parfor');  Opt.parfor  = 'off';  end
if ~isfield(Opt, 'gpu');     Opt.gpu     = 'off';  end

phys_modes = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2));

%% Phase 1: Build LHS matrices (once per run)
needs_beta = any(StabRes.betavec ~= 0);
if Opt.mmry == 1 && needs_beta && (isempty(MLA3) || isempty(MLA4))
    fprintf(' MLA3/MLA4 missing for 3D run — rebuilding LHS.\n');
    Opt.mmry = 0;
end
if Opt.mmry ~= 1
    Opt.mmry = 1;
    fprintf(' Constructing LHS matrices ...\n');
    tic
    [MLAMFD, MLA1, MLA2, MLA3, MLA4] = LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
    MFD_diff = MLAMFD - MLA1;
    clear MLAMFD
    fprintf(' LHS done in %.1f s\n', toc);
end

%% Phase 2: Pre-build LU factorisations (sequential — modifies Opt.LU_cache)
if ~isfield(Opt, 'LU_cache'); Opt.LU_cache = struct(); end

for j = phys_modes
    omega   = StabRes.omegavec(j);
    beta    = StabRes.betavec(j);
    ML_key  = ml_key(omega, beta);

    if ~isfield(Opt.LU_cache, ML_key)
        ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MFD_diff);

        if strcmpi(Opt.gpu, 'on')
            % GPU mode: store gpuArray(ML); \ refactorises on GPU each solve.
            % No decomposition caching — GPU does not support it for sparse complex.
            fprintf(' ML (%d,%d) [GPU] ...', Stab.Mvec(j), Stab.Nvec(j));
            Opt.LU_cache.(ML_key) = gpuArray(ML);
        else
            switch Opt.lu_mode
                case 'full'
                    fprintf(' LU (%d,%d) [cache all] ...', Stab.Mvec(j), Stab.Nvec(j));
                    Opt.LU_cache.(ML_key) = decomposition(ML, 'lu');

                case 'auto'
                    if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = 3; end
                    if ~isfield(Opt, 'lu_order');     Opt.lu_order = {}; Opt.lu_order_label = {}; end
                    n_cached = numel(fieldnames(Opt.LU_cache));
                    if n_cached >= Opt.lu_max_cache
                        old_key  = Opt.lu_order{1};
                        old_mode = Opt.lu_order_label{1};
                        Opt.LU_cache = rmfield(Opt.LU_cache, old_key);
                        Opt.lu_order(1) = [];
                        Opt.lu_order_label(1) = [];
                        fprintf(' LU (%d,%d) [replacing %s, %d cached] ...', Stab.Mvec(j), Stab.Nvec(j), old_mode, Opt.lu_max_cache);
                    else
                        fprintf(' LU (%d,%d) [%d/%d cached] ...', Stab.Mvec(j), Stab.Nvec(j), n_cached+1, Opt.lu_max_cache);
                    end
                    Opt.LU_cache.(ML_key) = decomposition(ML, 'lu');
                    Opt.lu_order{end+1}       = ML_key;
                    Opt.lu_order_label{end+1} = sprintf('(%d,%d)', Stab.Mvec(j), Stab.Nvec(j));

                case 'none'
                    fprintf(' LU (%d,%d) [direct] ...', Stab.Mvec(j), Stab.Nvec(j));
                    Opt.LU_cache.(ML_key) = ML;
            end
        end
        fprintf(' done\n');
        clear ML
    end
end

%% Phase 3: Solve — parfor or for over physical modes
nphys  = numel(phys_modes);
q_out  = zeros(nphys, size(StabRes.q, 2));

% Extract LU objects into a cell array (parfor requires indexing, not dynamic field names)
LU_cell = cell(nphys, 1);
for k = 1:nphys
    j = phys_modes(k);
    LU_cell{k} = Opt.LU_cache.(ml_key(StabRes.omegavec(j), StabRes.betavec(j)));
end

% Broadcast read-only copies for parfor
Stab_f   = Stab.f;    % rows indexed by mode j — each j independent
Mvec     = Stab.Mvec;
Nvec     = Stab.Nvec;

if strcmpi(Opt.parfor, 'on')
    fprintf(' Solving %d modes in parallel ...\n', nphys);
    tic_solve = tic;
    parfor k = 1:nphys
        j   = phys_modes(k);
        R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
        f_j = Stab_f(j,:);
        f_j(R_j ~= 0) = 0;
        q_out(k,:) = (LU_cell{k} \ (R_j + f_j.')).';   %#ok<PFBNS>
    end
    fprintf(' Solve done in %.1f s\n', toc(tic_solve));
elseif strcmpi(Opt.gpu, 'on')
    fprintf(' Solving %d modes on GPU ...\n', nphys);
    tic_solve = tic;
    for k = 1:nphys
        j   = phys_modes(k);
        R_j = RHS(Bound, StabGrid, Mode, Opt, StabRes, j);
        f_j = Stab_f(j,:);
        f_j(R_j ~= 0) = 0;
        rhs_gpu    = gpuArray(R_j + f_j.');
        q_out(k,:) = gather(LU_cell{k} \ rhs_gpu).';
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
        q_out(k,:) = (LU_cell{k} \ (R_j + f_j.')).';
        fprintf(' ... %.1fs\n', toc(tic_j));
    end
end

%% Collect results (including conjugate modes)
for k = 1:nphys
    j = phys_modes(k);
    StabRes.q(j,:)              = q_out(k,:);
    StabRes.q(Mode.nf+1-j,:)   = conj(q_out(k,:));
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

function ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MFD_diff)
if omega == 0 && beta == 0
    ML = MLA1 + MFD_diff;
elseif omega == 0
    ML = MLA1 + MLA3*beta + MLA4*beta^2;
elseif beta == 0
    ML = MLA1 + MLA2*omega;
else
    ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2;
end
end
