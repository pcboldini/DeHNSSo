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
% Stability grid for the DeHNSSo-v1 TS Blasius test case.
%
% Grid of the v1 caller
%   Callers/TS-Waves in a Blasius BL/Example_Caller_Blasius_BL_TS.m
% Grid.nx = 800 equidistant, Grid.ny = 40 Chebyshev, Grid.H = max(BF.Y),
% Grid.y_i = Grid.H/20, built from the v1 base flow "Data Files/BF_Blasius.mat".
% That base flow is not shipped with DeHNSSo-v2. Clone DeHNSSo-v1 next to this
% repository,
%   git clone https://github.com/SvenWesterbeek/DeHNSSo DeHNSSo-v1
% or set v1dir below. The resulting grid is shipped as
% sample/input/DeHNSSo_input_sample.mat.
%
% Writes input/DeHNSSo_input_v1_blasius.mat.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(fullfile(rootdir, 'gridgen'), genpath(fullfile(rootdir, 'gridgen', 'src')))

v1dir = fullfile(rootdir, '..', 'DeHNSSo-v1');   % DeHNSSo-v1 repository

%% Input
input.folder     = fullfile(v1dir, 'Data Files');
input.filename   = 'BF_Blasius.mat';
input.format     = 'mat';
input.structured = true;
if ~isfile(fullfile(input.folder, input.filename))
    error('gridgen_ts_blasius:noBF', 'DeHNSSo-v1 base flow not found: %s', fullfile(input.folder, input.filename));
end

%% Resampling — v1 values
bf  = load(fullfile(input.folder, input.filename));
BF  = bf.BF;
params.n_eta_new = 40;
params.n_xi_new  = 800;
params.H         = max(BF.Y(:));   % v1: Grid.H = max(BF.Y)
params.y_i       = params.H/20;    % v1: Grid.y_i = Grid.H/20
params.rescale   = false;
params.plot      = false;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, 'input');
output.filename = 'DeHNSSo_input_v1_blasius.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
