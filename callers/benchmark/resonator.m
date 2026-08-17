%% DeHNSSo — Helmholtz resonator in a compressible Blasius BL (TS-wave control)
%  A 1 mm Helmholtz-resonator neck is embedded in the wall of a Ma=0.061 Blasius
%  BL (Michelis 2023: U0 = 20.9 m/s, f_r = 295 Hz, x0 = 0.55 m). The resonator
%  imposes the reactive wall relation v'(x) - Y*p'(x) = 0 over the mouth (Opt.hr;
%  apply_resonator_wall). The TS wave is solved with a rigid wall (baseline) and
%  with the resonator, and their N-factor increment dN is compared with Michelis.
%
%  Base flow / grid: run baseflow/caller/caller_CBL_resonator.m, then
%  gridgen/benchmark/resonator.m (writes DeHNSSo_input_resonator.mat).

%% Setup
clear all; close all; clc
restoredefaultpath
rootdir = fullfile(fileparts(mfilename('fullpath')), '..', '..');
addpath(genpath(fullfile(rootdir, 'src')))
addpath(genpath(fullfile(rootdir, 'gridgen', 'src')))   % deriv2_uneven etc. used by init()

input.folder = fullfile(rootdir, 'input');
load(fullfile(input.folder, 'DeHNSSo_input_resonator.mat'));   % neck-clustered grid

Opt.compressible = 'on';        % (u,v,w,rho,T) path

% Michelis (2023) reference scales & neck geometry
U0   = 20.9;   lref = 2.6955e-4;   Re0 = 371;
f_TS = 295;                        % TS frequency = resonator resonance [Hz] (on resonance)
s_nd = 1e-3/lref;                  % neck (mouth) width [non-dim]
Xc   = 0.55/lref + s_nd/2;         % neck centre [non-dim]

%% Spectral content — single TS mode
Stab.M = 1;                        % omega (temporal) modes
Stab.N = 0;                        % beta  (spanwise) modes
Stab.omega_0 = 2*pi*f_TS*lref/U0;  % fundamental angular frequency (f_TS = 295 Hz)
Stab.beta_0  = 0;                  % 2D

%% Inflow — compressible local stability (slow Mack/TS mode)
Stab.IC         = 'CLST';
Stab.phaseRef   = 'pwall';
Opt.mode_select = 'slow';
Opt.linear      = 'on';            % linear: no NLT, single solve
Stab.A0_fund    = 1e-4;            % arbitrary in linear mode

%% Boundary conditions — no-slip isothermal wall (the neck row is set by Opt.hr)
Opt.bc_top          = {'H_DR', 'H_DR', 'H_DR', 'H_DR'};
Opt.bc_wall_thermal = 'isothermal';

%% Buffer — 'on' (amplitude damping) is cleaner than 'para' at this low Mach
Opt.buffer = 'on';
Opt.xb     = 90;                   % buffer start [% of domain]

%% Solver
Opt.plot    = 'off';               % suppress per-run plots (custom comparison below)
Opt.lu_mode = 'lapack_band';       % banded LU (compressible); build via mexbuild_dgbsv.m
Opt.parfor  = 'off';
Opt.solver  = 'picard';

%% Resonator — admittance imposed over the neck mouth
Ynd  = -12.7688;                   % non-dim resonator admittance at f_TS = 295 Hz (on resonance)
x    = StabGrid.xun(1, :);   nx = numel(x);
in_neck = (x >= Xc - s_nd/2) & (x <= Xc + s_nd/2);   % mouth stations
fprintf('Ynd = %.4f,  neck = %d stations\n', Ynd, nnz(in_neck));

%% Run — rigid-wall baseline, then the resonator (each a fresh operator)
Opt.rerun       = 'off';
Opt.save_output = [];

tic_total = tic;
StabRes0 = main(Stab, StabGrid, Opt);                          % baseline (rigid wall)

Opt.hr   = struct('Y', Ynd, 'mask', in_neck, 'x', x, 'taper', 0, 'uslip', false);
StabResA = main(Stab, StabGrid, Opt);                          % resonator
fprintf('\nTotal execution time: %.1f s\n', toc(tic_total));

%% Post-processing — N-factor increment dN (Michelis Eq. 4)
%  dN_u' at the u-peak, dN_v' at the v-peak, measured just upstream of the buffer.
jA    = find(abs(StabRes0.omegavec(:).' - Stab.omega_0) < 1e-9, 1);
imeas = max(1, round(Opt.xb*nx/100) - 1);
Rx    = sqrt(Re0 * x);

Au0 = StabRes0.A(jA, :);   AuA = StabResA.A(jA, :);
Av0 = Au0 .* max(abs(squeeze(StabRes0.v(jA,:,:))), [], 1);
AvA = AuA .* max(abs(squeeze(StabResA.v(jA,:,:))), [], 1);
dNu = log(AuA(imeas) / Au0(imeas));
dNv = log(AvA(imeas) / Av0(imeas));

fprintf('\n neck Rx = %.0f,  measure Rx = %.0f\n', sqrt(Re0*Xc), Rx(imeas));
fprintf(' dN_u'' = %+.3f,  dN_v'' = %+.3f    (Michelis ~ +0.4, amplifying)\n', dNu, dNv);

% N-factor, baseline vs resonator
xr  = 1:imeas;
Nu0 = log(Au0 ./ Au0(1));   NuA = log(AuA ./ AuA(1));
figure('Position', [100 100 760 460]);
plot(Rx(xr), Nu0(xr), 'k-',  'LineWidth', 1.5, 'DisplayName', 'baseline (rigid wall)'); hold on
plot(Rx(xr), NuA(xr), 'r--', 'LineWidth', 1.5, 'DisplayName', 'resonator');
xline(sqrt(Re0*Xc), 'b:', 'neck', 'HandleVisibility', 'off');
xlabel('$R_x = \sqrt{Re_0\,x}$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$N_{u''}$', 'Interpreter', 'latex', 'FontSize', 12)
title(sprintf('TS N-factor, Helmholtz resonator (dN_{u''} = %+.2f)', dNu))
legend('Location', 'best'); grid on; box on
