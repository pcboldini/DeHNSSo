%% DeHNSSo — Blasius boundary layer (TS waves)
% Tollmien-Schlichting wave growth in a zero-pressure-gradient flat-plate BL.
% Quick test: set nx = 800, ny = 40 in the input file (~1 h on a desktop PC).

%% Setup
clear all; close all; clc
restoredefaultpath
root = fileparts(mfilename('fullpath'));  % examples/ parent = repo root
addpath(genpath(fullfile(root, 'src')))

% Load stability grid from stability/ (produced by interpolation/main_interpol.m)
stability_folder = fullfile(root, 'stability');
load(fullfile(stability_folder, 'DeHNSSo_input.mat'));

% =========================================================================
%% Spectral content
% =========================================================================

Stab.N = 0;        % [-]  Spectral truncation: number of beta (spanwise) modes
Stab.M = 5;        % [-]  Spectral truncation: number of omega (temporal) modes
                   %      Total modes: (2N+1) x (2M+1), including conjugates

Stab.omega_0 = 0.034576;   % [-]  Fundamental angular frequency
                           %      (f = 90.58 Hz, F = 86.44 x 10^-6 at this Re)
Stab.beta_0  = 0;          % [-]  Fundamental spanwise wavenumber (0 = 2D)

% =========================================================================
%% Inflow conditions
% =========================================================================

Stab.IC       = 'ILST';   % Inflow method: 'ILST' (local stability) or 'LOAD'
Stab.phaseRef = 'pwall';  % Phase reference: 'pwall' (wall pressure) or 'umax'

% Initial amplitude vector — one entry per mode (size 1 x (2N+1)(2M+1)).
% Use ModeToModeNumber(m, n, Stab.M, Stab.N) to get the index for mode (m,n).
% Mode indexing example for M=1, N=1:
%      m:  -1   0   1
%   n= -1:  1   2   3
%   n=  0:  4   5   6
%   n=  1:  7   8   9
Stab.A0 = zeros(1, (2*Stab.N+1)*(2*Stab.M+1));
Stab.A0(ModeToModeNumber( 1, 0, Stab.M, Stab.N)) = 0.00125*sqrt(2);  % fundamental
Stab.A0(ModeToModeNumber(-1, 0, Stab.M, Stab.N)) = 0.00125*sqrt(2);  % conjugate

% =========================================================================
%% Boundary conditions — wall actuation and freestream forcing
% =========================================================================
% Prescribe inhomogeneous Dirichlet BCs at wall (y=0) and/or freestream (y=H).
% Both are commented out by default (homogeneous BCs apply).
%
% Stab.bcx         (1 x nx_bc)        x-locations of BC profile
% Stab.bcw(k,ix,j) (3 x nx_bc x nf)  wall BC:       k=1 u, k=2 v, k=3 w, mode j
% Stab.bct(k,ix,j) (3 x nx_bc x nf)  freestream BC: k=1 u, k=2 v, k=3 w, mode j
%
% Note: for the (0,0) MFD mode, the v-freestream row carries the transpiration
%       condition and is not overridden by bct — all other modes work normally.

% --- Wall actuation: Gaussian blowing/suction slot (TS wave control) ---
% j_fund = ModeToModeNumber( 1, 0, Stab.M, Stab.N);
% j_conj = ModeToModeNumber(-1, 0, Stab.M, Stab.N);
% x_dom  = StabGrid.xun;
% x0_act = x_dom(1) + 0.30*(x_dom(end) - x_dom(1));  % actuator centre
% sig    = 0.01*(x_dom(end) - x_dom(1));              % actuator half-width
% A_act  = 1e-4;                                       % actuation amplitude [-]
% Stab.bcx             = x_dom;
% Stab.bcw             = zeros(3, numel(x_dom), (2*Stab.N+1)*(2*Stab.M+1));
% Stab.bcw(2,:,j_fund) =  A_act * exp(-(x_dom - x0_act).^2 / sig^2);
% Stab.bcw(2,:,j_conj) = conj(Stab.bcw(2,:,j_fund));  % conjugate symmetry

% --- Freestream forcing: prescribe u at y=H (e.g. free-stream disturbance) ---
% j_fund = ModeToModeNumber( 1, 0, Stab.M, Stab.N);
% j_conj = ModeToModeNumber(-1, 0, Stab.M, Stab.N);
% x_dom  = StabGrid.xun;
% A_fst  = 1e-4;                                       % disturbance amplitude [-]
% Stab.bcx             = x_dom;
% Stab.bct             = zeros(3, numel(x_dom), (2*Stab.N+1)*(2*Stab.M+1));
% Stab.bct(1,:,j_fund) =  A_fst * ones(1, numel(x_dom));  % uniform u at freestream
% Stab.bct(1,:,j_conj) = conj(Stab.bct(1,:,j_fund));      % conjugate symmetry

% =========================================================================
%% Solver options
% =========================================================================

% --- Freestream BC ---
Opt.bc_top = 'dirichlet';  % 'dirichlet' u=v=w=0 at y=H (valid when domain is tall enough)
                            % 'neumann'   d/dy=0  at y=H (use when domain is truncated)

% --- Buffer ---
Opt.buffer = 'on';   % 'on'   amplitude damping (standard)
                     % 'para' parabolisation — zeroes streamwise diffusion (LPSE-like)
                     % 'off'  no buffer
Opt.xb     = 85;     % [-]  Buffer start as % of domain (default 85)

% --- Convergence ---
% Opt.Conv     = 1e-4;   % [-]  Convergence criterion (default 1e-4)
% Opt.ConvF    = 100;    % [-]  Relaxation factor during amplitude ramping (default 100)
% Opt.TH       = 1e-11;  % [-]  Nonlinear term activation threshold (default 1e-11)
% Opt.AFg      = 1.1;    % [-]  Amplitude factor growth rate (default 1.1)
% Opt.AMAX     = 0.1;    % [-]  Maximum amplitude factor (default 0.1)

% --- LU factorisation cache ---
Opt.lu_mode      = 'full';  % 'full' cache all LU factorisations (fastest, most memory)
                             % 'auto' cache up to lu_max_cache, evict oldest
                             % 'none' backslash each time (slowest, least memory)
% Opt.lu_max_cache = 3;      % max cached factorisations — only for 'auto' mode

% --- Parallel and GPU (both require Parallel Computing Toolbox) ---
Opt.parfor = 'off';    % 'on'  parfor over stations (LHS) and modes (solve) — auto-detected if omitted
Opt.gpu    = 'off';    % 'on'  GPU solve via gpuArray \ (beneficial for ny >> 80) — auto-detected if omitted

% --- Output ---
Opt.plot  = 'on';   % 'on'  intermediate amplitude/convergence plots during iteration
                    % 'off' suppress plots (useful for batch/benchmark runs)
% Opt.Sweep = 0;    % save intermediate results at each converged amplitude step

% =========================================================================
%% Run
% =========================================================================
% Set Opt.rerun = 'off' for a full run (builds LHS, returns it in LHS struct).
% Set Opt.rerun = 'on'  to reuse LHS from a previous run (skips assembly);
%   LU also reused when omega/beta are unchanged (actuation/amplitude sweeps),
%   recomputed when they change (frequency sweeps).
%
% Sweep example (frequency sweep):
%   omega_vec = Stab.omega_0 * linspace(0.8, 1.2, 5);
%   for k = 1:numel(omega_vec)
%       Stab.omega_0 = omega_vec(k);
%       if k == 1
%           [StabRes_sweep{k}, ~, LHS] = main(Stab, StabGrid, Opt);
%       else
%           Opt.rerun = 'on';
%           [StabRes_sweep{k}, ~, LHS] = main(Stab, StabGrid, Opt, LHS);
%       end
%   end

Opt.rerun = 'off';  % 'off' = full run | 'on' = reuse LHS from previous run

tic_total = tic;
if strcmpi(Opt.rerun, 'on')
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);
else
    [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
end
t_total = toc(tic_total);
w = whos;
fprintf('\n========================================\n');
fprintf('  Total execution time: %.2f s\n', t_total);
fprintf('  Workspace memory:     %.1f MB\n', sum([w.bytes])/1e6);
fprintf('========================================\n');

% =========================================================================
%% Post-processing
% =========================================================================
close all

% Amplitude plot — pass ref struct to overlay reference data, or omit it.
ref = struct('file', fullfile(stability_folder, 'literature', 'StabRes_Blasius_NPSE.mat'), 'label', 'NPSE');
plot_amplitudes(StabGrid, StabRes, ref);
% plot_amplitudes(StabGrid, StabRes);  % no reference
