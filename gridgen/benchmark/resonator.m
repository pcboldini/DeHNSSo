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
% Gridgen benchmark — Blasius boundary layer with a Helmholtz-resonator neck (Sec. 6.4)
%  Michelis et al. (2023) conditions: U0 = 20.9 m/s, resonance 295 Hz, a 1 mm
%  neck opening at x0 = 0.55 m. The 1 mm mouth has to carry enough stations for
%  the impedance row to be resolved, so the streamwise mesh is clustered on it
%  with a Gaussian point density: the target ratio is 4 at the neck, falling off
%  over 60 displacement thicknesses, which puts about 20 stations in the mouth.
%
%  Base flow: baseflow/benchmark/resonator.m (writes BF_compressible_M006_fine.mat
%  at high resolution, so this step is a clean downsampling).
%  Writes input/DeHNSSo_input_resonator.mat  (Table C.7: 5000 x 40, FD4)

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, ...
        genpath(fullfile(rootdir, 'src')), ...
        genpath(fullfile(rootdir, '..', 'src')))   % src/thermo/sutherland.m

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output', 'benchmark');
input.filename   = 'BF_compressible_M006_fine.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling
params.n_eta_new = 40;
params.n_xi_new  = 5000;
params.rescale   = false;
params.plot      = false;

%% Neck clustering — Gaussian point density centred on the mouth
% xi_nodes, like params.mug and params.spread, is in the units of the base flow's
% X field: displacement thicknesses here, so the 1 mm neck becomes 1e-3/lref.
BF   = load(fullfile(input.folder, input.filename)).BF;
lref = BF.lref;
s_nd = 1e-3/lref;                       % 1 mm neck, non-dimensional
Xc   = 0.55/lref + s_nd/2;              % neck centre
Xa   = BF.X(1);   Xb = BF.X(end);
rcl  = 4;                               % point-density ratio at the neck
wc   = 60;                              % Gaussian half-width [delta*]
Xr   = linspace(Xa, Xb, 120000);
dens = 1 + (rcl-1)*exp(-0.5*((Xr-Xc)/wc).^2);
Cc   = cumtrapz(Xr, dens);  Cc = (Cc-Cc(1))/(Cc(end)-Cc(1));
params.xi_nodes = interp1(Cc, Xr, linspace(0, 1, params.n_xi_new));
params.xi_nodes([1 end]) = [Xa Xb];

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

x = StabGrid.xun(1,:);
fprintf('\n%d x %d, %d stations in the 1 mm neck at R = %.0f\n', ...
        StabGrid.nx, StabGrid.ny, ...
        nnz(x >= Xc-s_nd/2 & x <= Xc+s_nd/2), sqrt(StabGrid.Re*Xc));
