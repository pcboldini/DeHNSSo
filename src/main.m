function [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS)
%MAIN  DeHNSSo main solver loop.
%
%   [StabRes, StabGrid]      = main(Stab, StabGrid, Opt)
%   [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt)       % also return LHS
%   [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS)  % reuse LHS (Opt.rerun='on')
%
%   Solves the nonlinear Harmonic Navier-Stokes equations iteratively:
%   1. Initialise grid, operators, inflow, and boundary conditions
%   2. Convergence loop: solve modes, update amplitudes, check convergence
%   3. Post-process: normalize shape functions, compute wavenumbers
%
%   LHS struct fields: MLA1, MLA2, MFD_diff — base LHS matrices that depend
%   only on the base flow, not on omega/beta. Pass back in with Opt.rerun='on'
%   to skip matrix assembly (useful for omega/beta/amplitude sweeps).

%% Input validation
if ~exist('StabGrid','var') || ~isstruct(StabGrid); StabGrid = struct(); end
if ~exist('Opt','var')      || ~isstruct(Opt);      Opt = struct();      end
if ~exist('Stab','var')     || ~isstruct(Stab);     Stab = struct();     end
if ~exist('LHS','var');     LHS = struct(); end
if ~isfield(Opt,'rerun');   Opt.rerun = 'off'; end

fprintf('\n');
fprintf(' ______     _   _  _   _  _____ _____       \n');
fprintf(' |  _  \\   | | | || \\ | |/  ___/  ___|      \n');
fprintf(' | | | |___| |_| ||  \\| |\\ `--.\\ `--.  ___  \n');
fprintf(' | | | / _ \\  _  || . ` | `--. \\`--. \\/ _ \\ \n');
fprintf(' | |/ /  __/ | | || |\\  |/\\__/ /\\__/ / (_) |\n');
fprintf(' |___/ \\___\\_| |_/\\_| \\_/\\____/\\____/ \\___/ \n');
fprintf('                                       v1.2  \n');
fprintf('  Delft Harmonic Navier-Stokes Solver\n');
fprintf('  =========================================\n');

%% Extract base-flow struct from StabGrid
fprintf('Extracting base flow ... ');
BF.Re  = StabGrid.Re;
BF.U   = StabGrid.U;    BF.V   = StabGrid.V;    BF.W   = StabGrid.W;
BF.dxU = StabGrid.dxU;  BF.dxV = StabGrid.dxV;  BF.dxW = StabGrid.dxW;
BF.dyU = StabGrid.dyU;  BF.dyV = StabGrid.dyV;  BF.dyW = StabGrid.dyW;
fprintf('done (9 fields, Re = %.2f)\n', BF.Re);

%% 1. Initialisation
[StabGrid, Init, Mode, Opt, Stab, StabRes] = init(StabGrid, Opt, Stab);
[Mode, StabRes] = inflow(Mode, Stab, StabGrid, StabRes, Opt);
[Bound] = boundary(StabGrid, Opt, Stab);

%% Print configuration (after init sets defaults)
print_config(StabGrid, Stab, Opt);

%% 2. Convergence loop
% Reuse LHS matrices from a previous run (Opt.rerun='on') or build fresh.
if strcmpi(Opt.rerun,'on') && isfield(LHS,'MLA1') && ~isempty(LHS.MLA1)
    MLA1     = LHS.MLA1;
    MLA2     = LHS.MLA2;
    if isfield(LHS, 'MLA3'); MLA3 = LHS.MLA3; else; MLA3 = []; end
    if isfield(LHS, 'MLA4'); MLA4 = LHS.MLA4; else; MLA4 = []; end
    MFD_diff = LHS.MFD_diff;
    Opt.mmry = 1;   % signal solve_modes to skip LHS construction
    % Restore cached LU factorisations so modes with unchanged omega/beta
    % (e.g. actuation sweep) skip refactorisation entirely.
    % For omega/beta sweeps the key changes and a new factorisation is computed.
    if isfield(LHS,'LU_cache') && ~isfield(Opt,'LU_cache')
        Opt.LU_cache = LHS.LU_cache;
    end
    fprintf(' Reusing LHS matrices from previous run.\n');
else
    MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MFD_diff = [];
end

fprintf('\n--- Convergence loop ---\n');
iter = 0;
t_loop = tic;

while max(abs(Opt.dal)) >= Opt.Conv || Opt.AF ~= 1
    iter = iter + 1;
    Opt.dal = zeros(Mode.nf, 1);

    % Solve for all active modes
    t_iter = tic;
    [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MFD_diff] = solve_modes( ...
        BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MFD_diff);

    % Damp inflow amplitude (nonlinear only)
    if nnz(Stab.A0) ~= 1 && Opt.k1 == 0
        [Opt, StabRes] = AVAL(StabGrid, Opt, Mode, StabRes);
    end

    % Update u,v,w,p, amplitudes, active modes, nonlinear forcing
    [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes);

    % Check convergence and ramp amplitude
    [Opt, Stab, StabRes] = check_convergence(StabGrid, Bound, Mode, Opt, Stab, StabRes);

    % Status (after check_convergence so dal is updated)
    n_active = numel(Mode.RunJ(Mode.RunJ >= round(Mode.nf/2)));
    max_A = max(StabRes.A(:));
    max_dal = max(abs(Opt.dal));
    fprintf('  Iter %3d | AF=%.3f | %d modes | max(A)=%.3e | max(dal)=%.3e | %.1f s\n', ...
        iter, Opt.AF, n_active, max_A, max_dal, toc(t_iter));

    % Intermediate plotting
    if strcmpi(Opt.plot,'on')
        plot_intermediate(StabGrid, StabRes, Mode, Stab, Opt);
    end

    % Sweep: save intermediate results
    if Opt.Sweep == 1 && max(abs(Opt.dal)) <= Opt.Sweep * Opt.Conv
        [Opt, StabRes] = save_sweep(StabGrid, Mode, Opt, Stab, StabRes);
    end

    % Convergence plot
    if strcmpi(Opt.plot,'on')
        plot_convergence(StabGrid, Mode, Stab, Opt);
    end

end
% Report peak memory (RSS)
try
    [~, rss] = system(sprintf('ps -o rss= -p %d', feature('getpid')));
    fprintf('--- Converged in %d iterations (%.1f s, peak %.1f GB) ---\n\n', ...
        iter, toc(t_loop), str2double(strtrim(rss))/1e6);
catch
    fprintf('--- Converged in %d iterations (%.1f s) ---\n\n', iter, toc(t_loop));
end

%% 3. Post-processing
StabRes = postprocess(StabGrid, Init, Mode, Opt, StabRes);

%% Pack LHS output (for reuse in subsequent runs)
LHS.MLA1     = MLA1;
LHS.MLA2     = MLA2;
LHS.MLA3     = MLA3;
LHS.MLA4     = MLA4;
LHS.MFD_diff = MFD_diff;
if isfield(Opt,'LU_cache')
    LHS.LU_cache = Opt.LU_cache;  % store factorisations for actuation/amplitude reruns
end

end
