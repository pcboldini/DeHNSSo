%% Gridgen benchmark — Compressible Blasius (FSC) base flow
%  Cartesian pass-through from baseflow/caller/caller_CBL.m. build_stab_grid
%  auto-detects rho/T, runs Sutherland, and forwards the thermo metadata.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, ...
        genpath(fullfile(rootdir, 'src')), ...
        genpath(fullfile(rootdir, '..', 'src')))   % src/thermo/sutherland.m

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output');
input.filename   = 'BF_compressible.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — pass-through
params.n_eta_new = [];
params.n_xi_new  = [];
params.rescale   = false;
params.plot      = false;

%% Grid operators (FD_xi_order_2 = 4: BF d2/dx2 viscous terms in LHS_compressible)
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 4;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_compressible.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
