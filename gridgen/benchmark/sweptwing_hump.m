%% Gridgen benchmark — SweptWing CFI over a smooth hump
%  Curvilinear structured input (80 x 1200, Y follows a Gaussian hump),
%  resampled onto Malik-Chebyshev eta x wall-arc-length xi.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_sweptwing_hump.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling
params.n_eta_new = 80;
params.n_xi_new  = 1200;
params.y_i       = [];
params.H         = [];
params.xi_range  = [];
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = false;
params.plot      = false;

%% Streamwise refinement — Gaussian cluster at the hump centre
x_m = 0.1837410879 / 2.1394e-4;  % hump centre ~859
params.mug    = x_m;
params.spread = 55.137096;       % half-width the refinement reaches (matches the heritage grid)
params.refine = 2;

%% Grid operators (fd eta: sparse pentadiagonal, avoids cheb OOM on laptop)
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_sweptwing_hump.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
