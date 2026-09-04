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
% Gridgen benchmark — Blasius BL base flow for the Newton benchmark
%  Resamples a Blasius base flow onto the 1050 x 60 grid of
%  callers/benchmark/blasius_newton.m: fourth-order finite differences in both
%  directions, R = 400-1200, H = 95.81.
%  Base flow: run baseflow/caller/caller_IBL.m with BL.L = 9*BL.S (domain end, R = 1200),
%  BL.H = 95.8144*lref, BL.y_i = BL.H/10, BL.nx = 2100, BL.ny = 101, saved as
%  baseflow/output/bf_blasius_newton.mat.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(rootdir, genpath(fullfile(rootdir, 'src')))

%% Input
input.folder     = fullfile(rootdir, '..', 'baseflow', 'output');
input.filename   = 'bf_blasius_newton.mat';
input.format     = 'mat';
input.structured = true;

%% Resampling
params.n_eta_new = 60;
params.n_xi_new  = 1050;
params.rescale   = false;
params.plot      = false;

%% Grid operators
params.FD_xi_order_1 = 4;
params.FD_xi_order_2 = 4;
params.FD_eta_method = 'fd';
params.FD_eta_order  = 4;

%% Output
output.folder   = fullfile(rootdir, '..', 'input');
output.filename = 'grid_blasius_fd_ny60_ext.mat';

%% Run
StabGrid = main_gridgen(input, params, output);
