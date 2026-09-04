%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Gridgen benchmark — SweptWing CFI over a smooth hump
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
