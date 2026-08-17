%% DeHNSSo — Synthetic jet actuator on a Blasius TS wave (benchmark)
%  Reference caller for the wall-actuator workflow. A linear TS wave (M=1) in a
%  Blasius BL is actuated by a wall-mounted synthetic jet: a wall-normal velocity
%  forcing v'|_wall over the orifice, at the TS frequency with prescribed
%  amplitude and phase. Models an active TS-control actuator (cf. the TS-control
%  study).
%
%  New actuator workflow (see docs/user_guide.md "Wall actuation"):
%    1. gridgen/benchmark/synthetic_jet.m clusters the mesh at the jet centre
%       (params.mug) and stores it as StabGrid.cluster.xc.
%    2. add_wall_forcing reads that centre and builds a raised-cosine forcing of a
%       PHYSICAL half-width w (not a grid-point count), so the refinement
%       resolves the same orifice better instead of shrinking it. Calls are
%       additive, so several actuators just stack.
%
%  Base flow / grid: run baseflow/.../bf_blasius.mat, then
%  gridgen/benchmark/synthetic_jet.m (writes DeHNSSo_input_jet.mat).

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()

input.folder = fullfile(rootdir, 'input');
load(fullfile(input.folder, 'DeHNSSo_input_jet.mat'));   % clustered grid + cluster.xc

%% Spectral content — single TS mode
Stab.M = 1;                 % omega (temporal) modes
Stab.N = 0;                 % beta  (spanwise) modes
Stab.omega_0 = 0.034576;    % fundamental angular frequency (F = 86e-6)
Stab.beta_0  = 0;           % 2D

%% Inflow — linear TS from local stability
Stab.IC       = 'ILST';
Stab.phaseRef = 'pwall';
Opt.linear    = 'on';       % linear: no NLT, single solve
Stab.A0_fund  = 0.0025*sqrt(2);

%% Freestream BC — free (unconstrained) v; u,w Dirichlet
Opt.bc_top = {'H_DR', 'free', 'H_DR'};

%% Buffer — 'on' (amplitude damping) | 'para' (parabolisation)
Opt.buffer = 'on';
Opt.xb     = 80;            % buffer start [% of domain]

%% Solver
Opt.plot        = 'on';
Opt.plot_metric = 'umax';
Opt.lu_mode     = 'full';
Opt.parfor      = 'off';
Opt.solver      = 'picard';

%% Baseline — incoming TS wave only (assembles and returns LHS for reuse)
Opt.rerun       = 'off';
Opt.save_output = [];
[StabRes0, StabGrid, LHS] = main(Stab, StabGrid, Opt);

%% Wall forcing: apply the profile gridgen stored (add_wall_forcing) ----------
% The shape, location and width are defined in gridgen/benchmark/synthetic_jet.m
% (StabGrid.forcing). Here we only set the strength. The conjugate mode is added
% automatically for nonlinear runs (here linear, so just the fundamental).
A_jet   = 2e-4;          % jet velocity amplitude (physical 2|u|, as Stab.A0)
phi_jet = 0;             % jet phase [rad] (control knob)

[Stab, jet] = add_wall_forcing(Stab, StabGrid, Opt, struct('A', A_jet, 'phi', phi_jet));
% To add a second forcing, call again — it superposes:
%   Stab = add_wall_forcing(Stab, StabGrid, Opt, struct('xc',x2,'w',12,'A',A2));

%% Actuated run — reuse LHS (only the RHS changes)
Opt.rerun = 'on';
[StabResA, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS);

%% Jet response at the orifice (mouth-averaged v and p at the wall)
vr = squeeze(StabResA.v(2,:,:)) .* StabResA.A(2,:);
pr = squeeze(StabResA.p(2,:,:)) .* StabResA.A(2,:);

% average over the orifice stations (jet.idx from add_wall_forcing)
pav  = mean(pr(end,jet.idx));
vav  = mean(vr(end,jet.idx));
vmax = max(vr(end,jet.idx));

ratav2 = abs(vav)/abs(pav)
ratav  = abs(vav/pav)
ratmax = abs(vmax)/abs(pav)
rata   = (angle(vav/pav))/pi
