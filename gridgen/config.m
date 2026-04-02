%% =========================================================
%% USER CONFIGURATION
%% =========================================================

%% Input file
input.format   = 'mat';               % 'mat', 'bin' or 'csv'
input.folder   = fullfile(fileparts(pwd), 'baseflow', 'output');
input.filename = '';       % specific file (leave empty = first file alphabetically)
input.z_slice    = 0;                   % z-value to extract from 3D OpenFOAM exports (leave empty = min(z))
input.structured = true;                % true: structured grid, false: unstructured scatter data

%% Malik-Chebyshev mesh preprocessing (create_ctype_mesh)
% Runs automatically for curvilinear meshes (detected from data).
% For Cartesian meshes, BF is used as-is.

% Grid dimensions (auto-detected; override if incorrect)
%   ctype_params.n_eta = 181;
%   ctype_params.n_xi  = 81;

% Target resolutions
ctype_params.n_eta_new = [64];        % wall-normal Chebyshev points
ctype_params.n_xi_new  = [137];        % streamwise  (leave empty = keep detected)

% Wall-normal domain height (leave empty = estimated from data)
ctype_params.H         = 1;

% Malik median collocation height (leave empty = H/3)
ctype_params.y_i       = 0.2;

% Trim streamwise domain as fraction of arc-length, e.g. [0.0, 0.9]
% (leave empty = full domain)
ctype_params.xi_range  = [];

% Non-dimensionalisation (required for dimensional input data)
% Leave empty if input is already non-dimensional (e.g. from BL solver).
ctype_params.Uref      = 20.9;        % reference velocity [m/s] (leave empty if already non-dim.)
ctype_params.lref      = 0.1; %6.0445e-4      % reference length   [m]  (leave empty if already non-dim.)
ctype_params.Re        = 400;       % Reynolds number Re = Uref * lref / nu

% Elliptic mesh smoothing (recommended for concave surfaces)
ctype_params.elliptic  = true;       % true: elliptic PDE smoothing, false: algebraic rays only
% ctype_params.elliptic_opts.max_iter   = 5000;
% ctype_params.elliptic_opts.tol        = 1e-6;
% ctype_params.elliptic_opts.omega      = 1.4;   % relaxation factor (1=none, >1=over-relax)
% ctype_params.elliptic_opts.n_uniform  = 201;   % uniform eta points for TDMA
% ctype_params.elliptic_opts.P_strength = 0;     % Poisson source strength (0 = Laplace)
% ctype_params.elliptic_opts.P_decay    = 5;     % source term decay rate

% Plots
ctype_params.plot            = true;
ctype_params.plot_xi         = 0;      % xi location for fig 5: fraction [0,1] or index (default: mid-domain)
ctype_params.plot_eta        = [0.9];  % eta location for fig 4: non-dim value (default: 5*delta_99)
ctype_params.plot_FD         = true;   % true: FD derivatives in fig 5, false: Chebyshev
ctype_params.plot_opts.xi_slices  = linspace(0, 1, 8);   % fractions [0,1] or indices
ctype_params.plot_opts.eta_slices = [0.0, 0.1, 0.3];     % fractions [0,1] or indices (leave empty to skip)

%% Grid generation options
StabGrid.Re    = 400;               % Reynolds number Re = Uref * lref / nu
StabGrid.plot  = true;

% Streamwise (xi) derivatives — always finite differences on the non-uniform xi grid
%   FD1_order : order of 1st derivative d/dxi  (used in metric tensors and base-flow xi-gradients)
%   FD2_order : order of 2nd derivative d2/dxi2 (used in metric tensors only; 2nd order is sufficient)
StabGrid.FD1_order = 4;             % 2, 4, or 6
StabGrid.FD2_order = 2;             % 2 or 4

% Wall-normal (eta) differentiation method — choose ONE:
%   'cheb' : spectral Chebyshev differentiation matrices (dense, O(N^2) LHS nonzeros)
%   'fd'   : sparse FD matrices on the Malik-Chebyshev grid (O(N) LHS nonzeros, ~80x fewer for N=100)
%            Same Malik-Chebyshev grid points in both cases; only D1/D2 differ.
%            FD_eta_order controls order for BOTH d/deta (D1) and d2/deta2 (D2).
StabGrid.eta_method   = 'cheb';       % 'cheb' or 'fd'
StabGrid.FD_eta_order = 4;          % FD order (only when eta_method='fd'): 2, 4, or 6

%% Output
output.folder   = fullfile(fileparts(pwd), 'stability');
output.filename = 'DeHNSSo_input.mat';
