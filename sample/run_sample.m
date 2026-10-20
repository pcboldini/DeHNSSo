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
% DeHNSSo sample run: TS waves in a Blasius boundary layer
%  Grid: 800 x 40, the grid of DeHNSSo-v1.
%  Reference: DeHNSSo-v1 solution on the same grid.

clear all; close all; clc
restoredefaultpath
sampledir = fileparts(mfilename('fullpath'));
rootdir   = fullfile(sampledir, '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()
load(fullfile(sampledir, 'input', 'DeHNSSo_input_sample.mat'));   % StabGrid

%% Spectral content — (2M+1)(2N+1) modes incl. conjugates
Stab.M = 5;                 % temporal modes
Stab.N = 0;                 % spanwise modes (0 = 2D)
Stab.omega_0 = 0.034576;    % fundamental frequency
Stab.beta_0  = 0;           % fundamental spanwise wavenumber

%% Inflow
Stab.IC       = 'ILST';     % 'ILST' (local stability) or 'LOAD'
Stab.phaseRef = 'pwall';    % 'pwall' (TS) or 'umax' (CFI)
Opt.linear    = 'off';      % 'on' linear (1 iter, NLT off) | 'off' nonlinear
Stab.A0_fund  = 0.0025*sqrt(2);   % fundamental amplitude (physical 2|u|)

%% Boundary conditions (freestream type per component u,v,w; wall is Dirichlet)
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};

%% Buffer — 'on' amplitude damping | 'para' parabolisation | 'off'
Opt.buffer = 'on';
Opt.xb     = 85;            % buffer start [% of domain]

%% Solver
Opt.plot        = 'off';
Opt.lu_mode     = 'full';   % 'full' | 'auto' | 'lapack_band' | 'none'
Opt.parfor      = 'off';
Opt.solver      = 'picard'; % 'picard' | 'newton'

%% Run
Opt.save_output = fullfile(rootdir, 'output', 'StabRes_sample.mat');   % [] to skip

tic_total = tic;
[StabRes, StabGrid] = main(Stab, StabGrid, Opt);
fprintf('\nTotal execution time: %.2f s\n', toc(tic_total));

%% Post-processing
load(fullfile(sampledir, 'reference', 'StabRes_sample_v1.mat'));   % StabResV1, StabGridV1
R   = sqrt(StabGrid.xun)   * sqrt(StabGrid.xun(1));
RV1 = sqrt(StabGridV1.xun) * sqrt(StabGridV1.xun(1));

figure('Position', [100 100 800 500])
h1 = semilogy(R, StabRes.A, 'k-', 'LineWidth', 1.2); hold on
h2 = semilogy(RV1(1:20:end), StabResV1.A(:,1:20:end), 'ro', 'MarkerSize', 5);
xlim(R([1 end])); ylim([1e-10 1e-1]); grid on
xlabel('$R=\sqrt{Re_x}$', 'Interpreter', 'latex')
ylabel('$\max|u''|$', 'Interpreter', 'latex')
legend([h1(1) h2(1)], {'DeHNSSo-v2', 'DeHNSSo-v1'}, 'Location', 'best')
