%% DeHNSSo — Blasius boundary layer (TS waves)
%  TS wave growth in a zero-pressure-gradient flat-plate BL.

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
Opt.linear  = 'on';        % 'on' linear (NLT off, 1 iter) | 'off' nonlinear
Stab.A0_fund = 0.00125*sqrt(2);

% =========================================================================
%% Boundary conditions
% =========================================================================
%
% Wall (y=0): always inhomogeneous Dirichlet.
%   Stab.bcw(k,ix,j) — wall BC value.  k=1 u, k=2 v, k=3 w.
%   Default: zero (no-slip).
%
% Freestream (y=H): per-component BC type.
%   Opt.bc_top{k} — type for component k=1 u, k=2 v, k=3 w.
%     'H_DR'  = homogeneous Dirichlet  (q = 0, default)
%     'IH_DR' = inhomogeneous Dirichlet (q = bct value, requires Stab.bct)
%     'H_NM'  = homogeneous Neumann    (dq/dy = 0)
%   Pressure: no BC at freestream (determined by the equations).
%
% Stab.bcx (1 x nx_bc) — x-locations on which bcw/bct are defined.
%   Profiles are interpolated onto the solver grid automatically.
%
% Both bcw and bct must be set on the same bcx. Always set conjugate modes.

Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};

% --- Example 1: wall blowing/suction slot ---
% j_fund = ModeToModeNumber( 1, 0, Stab.M, Stab.N);
% j_conj = ModeToModeNumber(-1, 0, Stab.M, Stab.N);
% x_dom  = StabGrid.xun;
% x0     = x_dom(1) + 0.30*(x_dom(end) - x_dom(1));
% sig    = 0.01*(x_dom(end) - x_dom(1));
% Stab.bcx             = x_dom;
% Stab.bcw             = zeros(3, numel(x_dom), (2*Stab.N+1)*(2*Stab.M+1));
% Stab.bcw(2,:,j_fund) = 1e-4 * exp(-(x_dom - x0).^2 / sig^2);
% Stab.bcw(2,:,j_conj) = conj(Stab.bcw(2,:,j_fund));

% --- Example 2: acoustic wave (v at freestream) ---
% Opt.bc_top           = {'H_DR', 'IH_DR', 'H_DR'};
% Stab.bcx             = x_dom;
% Stab.bct             = zeros(3, numel(x_dom), (2*Stab.N+1)*(2*Stab.M+1));
% Stab.bct(2,:,j_fund) = 1e-4 * ones(1, numel(x_dom));
% Stab.bct(2,:,j_conj) = conj(Stab.bct(2,:,j_fund));

% =========================================================================
%% Buffer — outflow treatment
% =========================================================================
%   'on'   = amplitude damping (scales solution -> 0 near outflow)
%   'para' = parabolisation (removes d2/dx2 in buffer -> no upstream reflections)
%   'off'  = no buffer
Opt.buffer  = 'para';
Opt.xb      = 85;             % buffer start [% of domain]

% =========================================================================
%% Solver options
% =========================================================================

Opt.plot    = 'on';         % 'on' | 'off'
% Opt.Conv  = 1e-4;         % convergence criterion
% Opt.TH    = 1e-11;        % NLT activation threshold

% =========================================================================
%% Performance
% =========================================================================

Opt.lu_mode = 'full';       % LU cache: 'full' | 'auto' (LRU) | 'none'
Opt.parfor  = 'off';        % parallel assembly/solve: 'on' | 'off'
Opt.gpu     = 'off';        % GPU sparse solve: 'on' | 'off'

% =========================================================================
%% Run
% =========================================================================

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

plot_stability(StabGrid, StabRes, Opt);
