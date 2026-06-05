%% Gridgen benchmark — Blasius BL base flow (cartesian CSV from OpenFOAM)

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_blasius.csv';
input.format     = 'csv';
input.structured = false;

%% Resampling parameters
params.n_eta_new = 60;                      % keep source ny
params.n_xi_new  = 1200;                      % keep source nx
params.rescale   = true;
params.Uref      = 20.9;            % CSV freestream U_e [m/s]
params.lref      = 2.891866e-04;   % Blasius length at the R=400 station (x=0.11569 m)
params.Re        = 400.03309;      % match bf_blasius.mat: R_inflow = sqrt(Re_x) = 400
params.plot      = true;

%% Streamwise window — reproduce the bf_blasius.mat R-window [400, 1067]
% The CSV spans x = [3.2e-4, 1.3] m, i.e. R = sqrt(Re_x) in [~21, ~1341] (leading
% edge -> far downstream). Keep only the sub-window R in [400, 1067], i.e.
% x in [0.1157, 0.8236] m, given as fractions of the x-span (flat plate: xi ~ x).
params.xi_range = [0.08877, 0.63344];
% params.xi_trim_inflow  = 0.02;     % drop first 2% of x
% params.xi_trim_outflow = 0.0;      % keep full outflow

% interp2 method used by resample_cartesian.
% Default 'makima' (handles non-uniform CFD spacing, no overshoot).
% 'spline' is smoother in the interior but can overshoot at boundaries.
% 'cubic' silently falls back to 'spline' on non-uniform grids.
% params.interp_method = 'makima';  % 'makima' | 'spline' | 'linear' | 'cubic'

%% Grid-generation options
params.FD_xi_order_1    = 4;
params.FD_xi_order_2    = 2;
params.FD_eta_method   = 'cheb';
params.FD_eta_order = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
