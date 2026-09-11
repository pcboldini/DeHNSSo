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
% DeHNSSo-v1 test case, reproduced with DeHNSSo-v2:
%   stationary crossflow instability over a forward-facing step (embedded boundary)
%
% Original v1 caller: Callers/CFI over a step in a swept-wing BL/Example_Caller_CFI_step_in_SweptWing_BL.m
% Reference data in this folder: StabRes_SweptWing_step_DNS.mat (DNS, paper coordinates)
%
% Run make_bf_cfi_step.m, then gridgen_cfi_step.m first.
% Linear, one stationary mode introduced by a blowing/suction strip upstream of the
% step (v1 Eq. 28). The comparison with the DNS is on the shape of u'_max(x): both are
% linear, so the HNS curve is scaled to the DNS at the crossflow peak.

clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))
load(fullfile(rootdir, 'input', 'DeHNSSo_input_v1_cfi_step.mat'));   % StabGrid

%% Spectral content
Stab.M       = 0;                               % stationary
Stab.N       = 1;                               % (0,0) and the +-beta pair
Stab.omega_0 = 0;
Stab.beta_0  = 2*pi * StabGrid.lref / 7.5e-3;   % lambda = 7.5 mm

%% Inflow: none. The mode is forced at the wall (below).
Stab.IC       = 'ILST';
Stab.phaseRef = 'umax';
Stab.A0       = zeros(1, (2*Stab.N+1)*(2*Stab.M+1));
Opt.linear    = 'on';

%% Wall actuator: blowing/suction strip on mode (0,1), window x = 135.5..143.6
nf    = (2*Stab.N+1)*(2*Stab.M+1);
j_act = ModeToModeNumber(0, 1, Stab.M, Stab.N);
xbcw  = [135.5, 143.6];
Abs   = 1.414*0.5e-5;
Stab.bcx = linspace(StabGrid.xun(1), StabGrid.xun(end), 5000);
Stab.bcw = zeros(3, numel(Stab.bcx), nf);
Stab.bct = zeros(3, numel(Stab.bcx), nf);
in_win = Stab.bcx >= xbcw(1) & Stab.bcx <= xbcw(2);
env    = 4*(Stab.bcx - xbcw(1)).*(xbcw(2) - Stab.bcx) / (xbcw(2)-xbcw(1))^2;
Stab.bcw(2, in_win, j_act) = Abs * env(in_win).^3;            % v at the wall

%% Boundary conditions and embedded boundary
Opt.bc_top = {'H_DR', 'H_DR', 'H_DR'};
Opt.solid    = 'auto';         % 'auto' masks on a snapped grid, ghost row otherwise | 'mask' | 'ghost'

%% Buffer
Opt.buffer   = 'para';
Opt.xb       = 85;
Opt.nltbufxb = 80;

%% Solver
Opt.plot        = 'on';
Opt.plot_metric = 'umax';
Opt.lu_mode     = 'auto';
Opt.lu_ordering = 'symrcm';  % 1500x150 stations: COLAMD fill grows with nx and needs
                             % ~32 GB for one factor; symrcm (banded) stays flat
Opt.parfor      = 'off';
Opt.solver      = 'picard';

%% Run
Opt.rerun       = 'off';
Opt.save_output = fullfile(fileparts(mfilename('fullpath')), 'output', 'StabRes_cfi_step.mat');
tic_total = tic;
[StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt);
fprintf('\nTotal execution time: %.2f s\n', toc(tic_total));

%% Post-processing: overlay on the DNS in paper coordinates
close all
x_paper = 3.604*StabGrid.xun + 220.6;
A_hns   = StabRes.A(end, :);                       % (0,1)
ref     = fullfile(fileparts(mfilename('fullpath')), 'StabRes_SweptWing_step_DNS.mat');
if isfile(ref)
    D = load(ref); X_dns = D.StabResDNS.X(:); A_dns = D.StabResDNS.A(:);
    [~, ip] = max(A_dns); scale = A_dns(ip) / interp1(x_paper, A_hns, X_dns(ip), 'linear', 'extrap');
else
    X_dns = []; A_dns = []; scale = 1;
end
if strcmpi(Opt.plot, 'on')
    figure(1); clf; plot(x_paper, scale*A_hns, 'k-', 'LineWidth', 1.5); hold on
    if ~isempty(X_dns); plot(X_dns, A_dns, 'm:', 'LineWidth', 1.5); legend({'DeHNSSo-v2 (ghost-cell EBM)','DNS'}); end
    xline(3.604*StabGrid.StepX + 220.6, 'k--'); grid on; xlim([700 920]); xlabel('x (paper)'); ylabel('u''_{max}');
end
