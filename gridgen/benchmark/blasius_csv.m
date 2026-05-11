%% Gridgen benchmark — Blasius BL base flow (cartesian CSV from OpenFOAM)

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'midPlane_0_1.csv';
input.format     = 'csv';
input.structured = false;

%% Resampling parameters
params.n_eta_new = [];                      % keep source ny
params.n_xi_new  = [];                      % keep source nx
params.rescale   = true;
params.Uref      = 20.9;
params.lref      = 2.688803650654011e-04;
params.Re        = 3.719126161394363e+02;
params.plot      = true;

%% Optional: streamwise trim & interpolation method (cartesian path)
% Drop fractions of the x-span at inflow / outflow before resampling.
% Useful for CFD sources with noisy inflow (BC adjustment) or outflow regions.
% params.xi_trim_inflow  = 0.02;     % drop first 2% of x
% params.xi_trim_outflow = 0.0;      % keep full outflow

% interp2 method used by resample_cartesian.
% Default 'makima' (handles non-uniform CFD spacing, no overshoot).
% 'spline' is smoother in the interior but can overshoot at boundaries.
% 'cubic' silently falls back to 'spline' on non-uniform grids.
% params.interp_method = 'makima';  % 'makima' | 'spline' | 'linear' | 'cubic'

%% Grid-generation options
params.FD1_order    = 4;
params.FD2_order    = 2;
params.eta_method   = 'cheb';
params.FD_eta_order = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
