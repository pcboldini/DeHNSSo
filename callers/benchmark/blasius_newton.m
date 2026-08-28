%% DeHNSSo — Blasius K-type (Scholten) via Newton pseudo-arclength continuation
%  Continuation in the inflow amplitude from A0_fund = 2.40e-3 to 5.657e-3, i.e.
%  u_rms = 0.004, on the FD4 ny=60 nx=1050 grid (R = 400-1200). A direct solve at
%  u_rms = 0.004 does not converge from a cold start, which is why this is a
%  continuation rather than a single Newton solve (cf. blasius.m, a fixed-amplitude
%  Picard run of a different case).
%
%  MEASURED on lsl, 27 Aug 2026, with the settings below:
%    206 steps, 44 factorisations, 461.6 min (7.7 h), peak ~80 GB resident.
%    Factor reuse is what makes it affordable: a step that factorises costs 7.0 min,
%    a step reusing the factor costs 0.1 min (~70x). Reuse also caps the arclength
%    step, though: ds only grows when the corrector takes <= newton_it_easy (4)
%    iterations, which happens on ~26% of steps, so ds settles near 0.8 and never
%    approaches newton_dsmax = 3. A smaller newton_reuse may finish sooner overall
%    despite more factorisations; untested.
%
%  VERIFIED against the M=4 run of 19 Jul 2026 at the same physical amplitude:
%    (0,0) 3.0%, (1,0) 3.6%, (2,0) 15.9%. The (2,0) gap is truncation, not error:
%    it is the highest mode M=2 keeps, so it absorbs what M=4 passes to (3,0)/(4,0).
%    By R ~ 1100 those two reach the same order as (2,0), so use M=4 if the top of
%    the domain matters.
%
%  M=4 variant: set Stab.M = 4. The Jacobian goes 1.512M -> 2.520M and the in-core
%  factor no longer fits (218 GB at relax=100, 328 GB at relax=200 on a 188 GB box),
%  so it also needs Opt.newton_mumps_ooc = 1. Budget ~19 h.
%
%  REQUIRES A WORKSTATION. At M=2 the dense front alone reaches ~80 GB, and
%  out-of-core does not help: MUMPS streams the stored factor, never the transient
%  front. Verified unrunnable on a 32 GB laptop at any M on this grid.

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
Stab.M = 2;                 % number of omega (temporal) modes
Stab.N = 0;                 % number of beta (spanwise) modes
Stab.omega_0 = 0.034576;    % fundamental angular frequency
Stab.beta_0  = 0;           % fundamental spanwise wavenumber (0 = 2D)

%% Inflow conditions
Stab.IC       = 'ILST';
Stab.phaseRef = 'pwall';
Opt.linear    = 'off';      % nonlinear
Stab.A0_fund  = 2.40e-3;    % continuation START amplitude

%% Boundary conditions / buffer
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};
Opt.buffer = 'para';        % parabolisation: no upstream reflections
Opt.xb     = 85;            % buffer start [% of domain]

%% Solver options
Opt.plot    = 'off';
Opt.rerun   = 'off';
Opt.solver  = 'newton';

%% Pseudo-arclength continuation (src/newton/solve_newton.m)
% A0_fund is the PHYSICAL peak 2|u_hat|, halved internally at the seed, so the
% u_rms = 0.004 target is 0.004*sqrt(2) = 5.657e-3. (Runs before the amplitude
% convention was unified used half this, 2.83e-3, for the same physical state.)
Opt.newton_target    = 5.657e-3;   % target A0_fund  (u_rms = 0.004)
Opt.newton_start     = 2.40e-3;
Opt.newton_seed2     = 2.50e-3;    % second seed point
Opt.newton_reuse     = 8;          % reuse the factor across steps; see header

% Ordering MUST follow the wall-normal discretisation -- each is fatal on the
% other operator (measured Jul 2026, lsl, filled M=4 Jacobian):
%   Chebyshev (dense-in-y cliques): AMD 123 GB OK  / METIS 216 GB -13
%   FD (sparse regular stencil)   : METIS 181 GB OK / AMD   311 GB -13
% 7 = auto; on the FD grid it selects METIS (verified identical to ICNTL(7)=5).
Opt.newton_mumps_ordering = 7;

% M=2 fits in core on a 188 GB box (~80 GB peak), so OOC is left off: it costs
% ~1.7x and buys nothing here. Enable it for M=4 (see header). If you do, point
% MUMPS_OOC_TMPDIR at a filesystem with room -- solve_newton refuses to start an
% out-of-core factorisation with less than 50 GB free.
% Opt.newton_mumps_ooc = 1;

% newton_it_thresh, newton_maxit, newton_dsmax and newton_tol are left unset:
% the solver defaults (8, 25, 3, 1e-6) are what this case used.

tag = sprintf('M%d_fd_ny60', Stab.M);
Opt.newton_seedfile = fullfile(rootdir, 'output', sprintf('arclength_%s_seed.mat', tag));
Opt.newton_ckfile   = fullfile(rootdir, 'output', sprintf('arclength_%s_state.mat', tag));
Opt.newton_resume   = 'off';       % 'on' resumes from newton_ckfile

%% Run
tic_total = tic;
[StabRes, info] = solve_newton(Stab, StabGrid, Opt);
t_total = toc(tic_total);

fprintf('\n========================================\n');
fprintf('  status:       %s\n', info.status);
fprintf('  A0 reached:   %.6e  (u_rms = %.5f)\n', info.lambda, info.lambda/sqrt(2));
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
