%% DeHNSSo — Blasius boundary layer (TS waves)
% TS wave growth in a zero-pressure-gradient flat-plate BL.

%% Setup
clear all; close all; clc
restoredefaultpath
root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(root, '..', 'src')))

stability_folder = fullfile(root, '..', 'stability');
load(fullfile(stability_folder, 'DeHNSSo_input.mat'));

% =========================================================================
%% Spectral content
% =========================================================================

Stab.M = 1;                % number of omega (temporal) modes
Stab.N = 0;                % number of beta (spanwise) modes
                            % total modes: (2M+1)(2N+1), including conjugates
Stab.omega_0 = 0.034576;   % fundamental angular frequency
Stab.beta_0  = 0;          % fundamental spanwise wavenumber (0 = 2D)

% =========================================================================
%% Inflow conditions
% =========================================================================

Stab.IC       = 'ILST';    % 'ILST' (local stability) or 'LOAD' (from file)
Stab.phaseRef = 'pwall';   % phase reference: 'pwall' (wall pressure) or 'umax'
% Stab.ICfile = fullfile(stability_folder, 'StabRes_previous.mat');  % required for LOAD

% Fundamental disturbance amplitude.
% The full amplitude vector (fundamental + conjugate) is built in init.
% For per-mode control, set Stab.A0 directly (size 1 x (2N+1)(2M+1)).
Opt.linear  = 'on';        % 'on' linear (NLT off, 1 iter) | 'off' nonlinear
Stab.A0_fund = 0.00125*sqrt(2);

% =========================================================================
%% Boundary conditions — wall actuation and freestream forcing
% =========================================================================
% Inhomogeneous Dirichlet BCs at wall (y=0) and/or freestream (y=H).
% Stab.bcx         (1 x nx_bc)        x-locations of BC profile
% Stab.bcw(k,ix,j) (3 x nx_bc x nf)  wall BC:       k=1 u, k=2 v, k=3 w
% Stab.bct(k,ix,j) (3 x nx_bc x nf)  freestream BC: k=1 u, k=2 v, k=3 w

% --- Wall actuation: Gaussian blowing/suction slot ---
% j_fund = ModeToModeNumber( 1, 0, Stab.M, Stab.N);
% j_conj = ModeToModeNumber(-1, 0, Stab.M, Stab.N);
% x_dom  = StabGrid.xun;
% x0_act = x_dom(1) + 0.30*(x_dom(end) - x_dom(1));
% sig    = 0.01*(x_dom(end) - x_dom(1));
% A_act  = 1e-4;
% Stab.bcx             = x_dom;
% Stab.bcw             = zeros(3, numel(x_dom), (2*Stab.N+1)*(2*Stab.M+1));
% Stab.bcw(2,:,j_fund) = A_act * exp(-(x_dom - x0_act).^2 / sig^2);
% Stab.bcw(2,:,j_conj) = conj(Stab.bcw(2,:,j_fund));

% --- Freestream forcing: prescribe u at y=H ---
% Stab.bct             = zeros(3, numel(x_dom), (2*Stab.N+1)*(2*Stab.M+1));
% Stab.bct(1,:,j_fund) = 1e-4 * ones(1, numel(x_dom));
% Stab.bct(1,:,j_conj) = conj(Stab.bct(1,:,j_fund));

% =========================================================================
%% Solver options
% =========================================================================

Opt.bc_top  = 'dirichlet';  % 'dirichlet' u=v=w=0 at y=H | 'neumann' d/dy=0
Opt.buffer  = 'para';         % 'on' amplitude damping | 'para' parabolisation | 'off'
Opt.xb      = 85;           % buffer start [% of domain]
Opt.lu_mode = 'full';       % 'full' | 'auto' | 'none'
Opt.parfor  = 'off';        % 'on' | 'off'
Opt.gpu     = 'off';        % 'on' | 'off'
Opt.plot    = 'on';         % 'on' | 'off'
% Opt.Conv  = 1e-4;         % convergence criterion
% Opt.TH    = 1e-11;        % NLT activation threshold

% =========================================================================
%% Run
% =========================================================================
% Opt.rerun = 'on' reuses LHS from a previous run (skips matrix assembly).

Opt.rerun = 'off';

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

if strcmpi(Opt.linear, 'on')
    plot_amplitudes(StabGrid, StabRes);
else
    ref = struct('file', fullfile(stability_folder, 'literature', 'StabRes_Blasius_NPSE.mat'), 'label', 'NPSE');
    plot_amplitudes(StabGrid, StabRes, ref);
end

% N-factor, growth rates, wavenumber, eigenfunctions
plot_stability(StabGrid, StabRes, Opt);

