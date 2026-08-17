%% Generate the compressible Blasius base flow for the resonator benchmark (Ma=0.061)
%  Michelis (2023) conditions: U0 = 20.9 m/s, adiabatic flat plate, R in [371, 3339],
%  H = 107. Written at high resolution (6000 x 101) so resampling onto the stability
%  grid in gridgen/benchmark/resonator.m is a clean downsampling. Same similarity
%  solver and continuity closure as caller_CBL.m. Writes
%  baseflow/output/BF_compressible_M006_fine.mat.

clear; clc
rootdir = fullfile(fileparts(mfilename('fullpath')), '..');
addpath(genpath(rootdir));

%% M006 parameters (Michelis 2023)
BL.Ma    = 0.0608974;  BL.Pr = 0.72;  BL.gamma = 1.4;
BL.Sstar = 0.376599;                       % 110.4/293 (Sutherland)
BL.beta_H = 0;  BL.We_ratio = 0;
BL.wall_bc = 'adiabatic';  BL.Tw = 1.0;
BL.eta_inf = 100;  BL.npts = 501;

V_inf   = 20.9038;   BL.nu = 1.51875e-5;
R_start = 371;   Xend = 3339;              % BF.X = x/lref runs R_start .. Xend
R_end   = sqrt(Xend * R_start);            % Reynolds R at outflow
BL.S    = R_start^2 * BL.nu / V_inf;
BL.L    = R_end^2   * BL.nu / V_inf - BL.S;
Uref    = V_inf;
lref    = sqrt(BL.S * BL.nu / Uref);
Re_0    = Uref * lref / BL.nu;
BL.H    = 107 * lref;   BL.y_i = 3.8 * lref;

BL.nx   = 6000;   BL.ny = 101;             % fine (clean downsampling in gridgen)
BL.Ue   = ones(1, BL.nx) * V_inf;   BL.We = 0;

%% Similarity solve + fields on the Malik-Chebyshev grid (mirror caller_CBL)
[~, BL] = compressibleFSC(BL);
nx = BL.nx;  ny = BL.ny;
x_phys = linspace(BL.S, BL.S + BL.L, nx);
[yvec, DM] = chebdif(ny, 2);
[y_phys, ~, ~] = MappingMalik(BL.H, BL.y_i, yvec', DM(:,:,1), DM(:,:,2));
eta_sim = BL.eta;

u = zeros(ny,nx); w = zeros(ny,nx); rho = zeros(ny,nx); T = zeros(ny,nx);
for i = 1:nx
    eta_local = y_phys * sqrt(BL.Ue(i) / (BL.nu * x_phys(i)));
    u(:,i)   = BL.Ue(i) * interp1(eta_sim, BL.fp,        eta_local, 'spline', BL.fp(end));
    w(:,i)   = BL.We    * interp1(eta_sim, BL.g,         eta_local, 'spline', BL.g(end));
    rho(:,i) = interp1(eta_sim, BL.rho_prof, eta_local, 'spline', BL.rho_prof(end));
    T(:,i)   = interp1(eta_sim, BL.theta,    eta_local, 'spline', BL.theta(end));
end
v = continuity_V(y_phys, x_phys, u, rho);      % compressible continuity (canonical)

%% Assemble non-dim BF struct
BF.Ma=BL.Ma; BF.Pr=BL.Pr; BF.gamma=BL.gamma; BF.Sstar=BL.Sstar;
BF.Uref=Uref; BF.lref=lref; BF.nu=BL.nu; BF.Re=Re_0;
BF.X = x_phys/lref;   BF.Y = (y_phys/lref).';
BF.H = BL.H/lref;     BF.y_i = BL.y_i/lref;
BF.U = u/Uref;  BF.V = v/Uref;  BF.W = w/Uref;  BF.rho = rho;  BF.T = T;

outfile = fullfile(rootdir, 'output', 'BF_compressible_M006_fine.mat');
save(outfile, 'BF', '-v7.3');
fprintf('wrote %s  (nx=%d, ny=%d, Re_0=%.1f, Ma=%.4f)\n', outfile, nx, ny, Re_0, BL.Ma);
fprintf('X=[%.0f %.0f], Y=[%.2f %.1f], H=%.1f ; max|T-1|=%.2e max|rho-1|=%.2e\n', ...
        BF.X(1), BF.X(end), min(BF.Y), max(BF.Y), BF.H, max(abs(BF.T(:)-1)), max(abs(BF.rho(:)-1)));
