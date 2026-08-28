%% Gridgen benchmark — flat-plate SweptWing CFI base flow
%  Curvilinear structured input (80 x 1200 body-fitted, Y constant per column),
%  resampled onto Malik-Chebyshev eta x wall-arc-length xi.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_sweptwing_flat.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling
params.n_eta_new = 40;
params.n_xi_new  = 1200;
params.y_i       = [];
params.H         = [];
params.xi_range  = [];
params.xi_trim_inflow  = [];
params.xi_trim_outflow = [];
params.rescale   = false;
params.plot      = false;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_sweptwing_flat.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
