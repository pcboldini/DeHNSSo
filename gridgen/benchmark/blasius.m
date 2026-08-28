%% Gridgen benchmark — Blasius BL base flow
%  Cartesian structured input, already at DeHNSSo resolution (pass-through).

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_blasius.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — pass-through (empty counts keep the source grid)
params.n_eta_new = [];
params.n_xi_new  = [];
params.rescale   = false;
params.plot      = false;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_blasius.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
