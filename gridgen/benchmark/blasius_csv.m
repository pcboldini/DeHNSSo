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
% Gridgen benchmark — Blasius BL base flow (cartesian CSV from OpenFOAM)

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'bf_blasius.csv';
input.format     = 'csv';
input.structured = false;

% base flow is stored gzipped in the repo; unpack on first use
if ~isfile(fullfile(input.folder, input.filename))
    gunzip(fullfile(input.folder, [input.filename '.gz']));
end

%% Resampling
params.n_eta_new = 60;
params.n_xi_new  = 1200;
params.rescale   = true;
params.Uref      = 20.9;
params.lref      = 2.891866e-04;   % Blasius length at R=400 (matches bf_blasius.mat)
params.Re        = 400.03309;
params.plot      = true;
params.xi_range  = [0.08877, 0.63344];   % R window [400, 1067] as x-span fractions

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 2;
params.FD_eta_method = 'cheb';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'DeHNSSo_input_blasius_csv.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
