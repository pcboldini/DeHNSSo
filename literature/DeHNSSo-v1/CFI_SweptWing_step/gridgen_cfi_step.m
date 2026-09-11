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
% Stability grid for the DeHNSSo-v1 forward-facing-step test case (embedded boundary).
%
% The wall stays flat; the step is detected from the base flow by gridgen and
% enters the solver as masks (setup_boundary).
%
% Base flow: bf_cfi_step_dns.mat, the spanwise mean of the three-dimensional
% DNS reduced by make_bf_from_dns.m to 2109 x 305 points over x = 125..220.
% Window x = 129.5..215.1 is the v1 paper Fig. 13 range (paper x = 3.604 x + 220.6,
% i.e. 685..993 in delta_0), step at x = 177.615, height 0.9656 from the DNS metadata.
%
% Uniform xi and a Malik wall-normal distribution with the median collocation height
% at the step top, which puts 38 of the 80 points below the step: the same step
% resolution the earlier 150-point grid reached, at a quarter of the factor memory.
%
% Writes input/DeHNSSo_input_v1_cfi_step.mat.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(fullfile(rootdir, 'gridgen'), genpath(fullfile(rootdir, 'gridgen', 'src')))

%% Input (trimmed base flow beside this file)
input.folder     = fileparts(mfilename('fullpath'));
input.filename   = 'bf_cfi_step_dns.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling
bf = load(fullfile(input.folder, input.filename));  BFr = bf.BF;
params.xi_range  = ([129.5, 215.1] - BFr.X(1)) / (BFr.X(end) - BFr.X(1));
params.n_xi_new  = 1500;
params.n_eta_new = 80;
params.y_i       = 1.0;    % median collocation height at the step top
params.H         = [];
params.rescale   = false;
params.plot      = false;

% Step geometry from the DNS metadata; detection confirms it independently.
params.StepType  = 'FFS';
params.StepX     = BFr.StepX;
params.StepH     = BFr.StepH;

%% Forward-facing step: nothing to set. gridgen reads the solid off the base flow
% (zero-velocity nodes above the wall), classifies it, snaps a station onto the
% face and a node onto the top, and stores it as StabGrid.step for setup_boundary.

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 4;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, 'input');
output.filename = 'DeHNSSo_input_v1_cfi_step.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
