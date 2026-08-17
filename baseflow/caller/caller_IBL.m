clc; clear all; format long
caller_dir = fileparts(mfilename('fullpath'));
rootdir    = fullfile(caller_dir, '..', '..');
addpath(genpath(rootdir));

BL.model = 0;                % 0 = incompressible, 1 = compressible

%% Inputs

% Mesh
BL.S   = 0.2418;             % streamwise domain start [m]
BL.L   = 1.7213;             % streamwise domain end [m]
BL.H   = 0.0585;             % wall-normal domain height [m]
BL.y_i = BL.H/20;            % Malik median node height [m]
BL.nx  = 800;                % streamwise gridpoints
BL.ny  = 40;                 % wall-normal gridpoints

% Freestream
V      = 10;                 % edge velocity at inlet [m/s]
BL.Ue  = ones(1,BL.nx)*V;    % streamwise edge velocity (ZPG: constant)
BL.We  = 0;                  % spanwise edge velocity (unswept: 0)
BL.nu  = 1.511e-5;           % kinematic viscosity [m^2/s]

plt.mainFig = 1;             % plot after solving (0 = off)

%% Boundary-layer fields
[BL] = incompressibleBL_B1(BL);

%% Reference scales and non-dimensionalisation
BF.Uref = V;                          % reference velocity [m/s]
BF.lref = sqrt(BL.S*BL.nu/BF.Uref);   % Blasius length scale [m]
BF.nu   = BL.nu;
BF.Re   = BF.Uref*BF.lref/BF.nu;      % Re_0 = sqrt(Re_x)

BF.X = BL.x/BF.lref;
BF.Y = BL.y/BF.lref;
BF.U = BL.u/BF.Uref;
BF.V = BL.v/BF.Uref;
BF.W = BL.w/BF.Uref;
% Derivatives are produced by gridgen (compute_bf_derivatives).

%% Postprocessing
[BL] = post_BL(BL);

%% Save
outdir = fullfile(caller_dir, '..', 'output');
if ~exist(outdir, 'dir'); mkdir(outdir); end
save(fullfile(outdir, 'BF_Blasius.mat'), 'BF');
