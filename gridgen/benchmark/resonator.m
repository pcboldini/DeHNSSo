%% Gridgen benchmark — compressible Blasius BL, mesh clustered at the resonator neck
%  Michelis (2023) conditions: U0 = 20.9 m/s, f_r = 295 Hz, 1 mm neck at x0 = 0.55 m.
%  The streamwise mesh is refined around the neck the same way it is around a hump
%  or a wall actuator: params.forcing marks the neck centre and gridgen clusters
%  the stations there so the 1 mm mouth carries enough nodes for the impedance BC.
%
%  Base flow: run baseflow/caller/caller_CBL_resonator.m first.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, ...
        genpath(fullfile(rootdir, 'src')), ...
        genpath(fullfile(rootdir, '..', 'src')))   % src/thermo/sutherland.m

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output');
input.filename   = 'BF_compressible_M006_fine.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling — downsample the fine source (clustering redistributes the stations)
params.n_eta_new = 40;
params.n_xi_new  = 4000;
params.rescale   = false;
params.plot      = false;

%% Neck clustering — same path as the hump / wall-forcing cases.
% params.forcing marks the neck centre; gridgen refines the streamwise mesh around
% it so the 1 mm mouth carries enough nodes for the impedance BC. StabGrid.forcing
% is produced but unused (the resonator drives the wall through Opt.hr, not
% add_wall_forcing).
lref = 2.6955e-4;   s_nd = 1e-3/lref;   Xc = 0.55/lref + s_nd/2;   % Michelis neck [non-dim]
params.forcing = struct('xc', Xc, 'w', 30, 'shape', 'cosine');
params.ag      = 0.2;       % peak/edge spacing ratio (smaller = tighter neck cluster)

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 4;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_resonator.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
