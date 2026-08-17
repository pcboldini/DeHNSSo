%% DeHNSSo — Blasius K-type (Scholten) via Newton pseudo-arclength continuation
%  Reproduces the lsl production run: continuation in the inflow amplitude
%  A0_fund from 2.40e-3 up to 5.657e-3 (u_rms = 0.004), using the Newton
%  arclength solver with cross-step MUMPS factor reuse (src/newton/solve_newton.m).
%
%  A direct solve at u_rms=0.004 does not converge from a cold start, which is
%  why this case is a continuation rather than a single Newton solve (cf.
%  blasius.m, which is a fixed-amplitude Picard run of a different case).

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()

input.folder = fullfile(rootdir, 'input');
gridname     = 'grid_blasius_fd_ny60_ext.mat';   % FD4 ny=60, nx=1050, R=400-1200
G = load(fullfile(input.folder, gridname)); StabGrid = G.StabGrid;

%% Spectral content
%  MEMORY: the coupled Jacobian is factorised whole. The lsl runs peaked near
%  ~120 GB RSS at M=3-4 on this grid, so M>=3 needs the workstation. Reduce M
%  (or the grid) to fit a laptop.

Stab.M = 4;                 % number of omega (temporal) modes
Stab.N = 0;                 % number of beta (spanwise) modes
Stab.omega_0 = 0.034576;    % fundamental angular frequency
Stab.beta_0  = 0;           % fundamental spanwise wavenumber (0 = 2D)

%% Inflow conditions

Stab.IC       = 'ILST';
Stab.phaseRef = 'pwall';
Opt.linear    = 'off';      % nonlinear
Stab.A0_fund  = 2.40e-3;    % continuation START amplitude (lambda_0)

%% Boundary conditions / buffer

Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};
Opt.buffer = 'para';        % parabolisation: no upstream reflections
Opt.xb     = 85;            % buffer start [% of domain]

%% Solver options

Opt.plot    = 'off';
Opt.lu_mode = 'auto';
Opt.parfor  = 'off';
Opt.rerun   = 'off';

Opt.solver              = 'newton';
Opt.newton_precond      = 'block';
Opt.newton_block_factor = 'lu';

%% Pseudo-arclength continuation (src/newton/solve_newton.m)
% Corrector holds ONE coupled MUMPS factor across steps (back-solves only) and
% refactors on age/stall. That reuse is the speedup: on lsl this reached the
% target in 3.2 h / 69 factors, vs 10.1 h / 100 factors refactoring every step.

Opt.newton_target    = 5.657e-3;   % target A0_fund  (u_rms = 0.004)
Opt.newton_start     = 2.40e-3;    % continuation start
Opt.newton_seed2     = 2.50e-3;    % second seed point
Opt.newton_reuse     = 0;          % refactor every step via the safe path

% Ordering MUST follow the wall-normal discretisation -- each is fatal on the
% other operator (measured Jul 2026, lsl, filled M=4 Jacobian):
%   Chebyshev (dense-in-y cliques): AMD 123 GB OK  / METIS 216 GB -13
%   FD (sparse regular stencil)   : METIS 181 GB OK / AMD   311 GB -13
% 7 = auto; on the FD grid it selects METIS (verified identical to ICNTL(7)=5).
Opt.newton_mumps_ordering = 7;

% This grid does not fit in core: measured Jul 2026 the filled Jacobian asked
% 218 GB (relax=100) then 328 GB (relax=200) on a 188 GB box. Skip straight to
% out-of-core rather than burning a failed analysis + pause(20) per attempt on
% every one of ~120 steps. Spills to /tmp (842 GB free); costs ~1.7x.
Opt.mumps_ooc = 1;
Opt.newton_it_thresh = 8;          % refactor next step if corrector exceeded this
Opt.newton_maxit     = 25;         % max corrector back-solves before stale
Opt.newton_dsmax     = 3;          % max arclength step
Opt.newton_tol       = 1e-6;       % corrector residual tolerance

tag = sprintf('M%d_fd_ny60', Stab.M);
Opt.newton_seedfile = fullfile(rootdir, 'output', sprintf('arclength_%s_seed.mat', tag));
Opt.newton_ckfile   = fullfile(rootdir, 'output', sprintf('arclength_%s_state.mat', tag));
Opt.newton_resume   = 'off';       % from scratch: no checkpoint is valid on this grid

%% Run

tic_total = tic;
[StabRes, info] = solve_newton(Stab, StabGrid, Opt);
t_total = toc(tic_total);

fprintf('\n========================================\n');
fprintf('  status:       %s\n', info.status);
fprintf('  A0 reached:   %.6e  (u_rms = %.5f)\n', info.lambda, sqrt(2)*info.lambda);
fprintf('  steps:        %d\n', info.steps);
fprintf('  factorisations: %d\n', info.nfac_total);
fprintf('  wall time:    %.1f min (%.2f h)\n', t_total/60, t_total/3600);
fprintf('========================================\n');

save_output = fullfile(rootdir, 'output', sprintf('StabRes_blasius_newton_%s.mat', tag));
save(save_output, 'StabRes', 'StabGrid', 'info', '-v7.3');
fprintf('  saved: %s\n', save_output);

%% Post-processing
close all
if strcmpi(Opt.plot, 'on')
    plot_amplitudes(StabGrid, StabRes, [], [], 'urms');
    plot_stability(StabGrid, StabRes, Opt);
end
