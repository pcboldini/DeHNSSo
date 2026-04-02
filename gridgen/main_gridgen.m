%% =========================================================
%  Interpolation pipeline for DeHNSSo
%  Edit config.m to set all parameters.
%% =========================================================

clear all;  close all;  clc;
restoredefaultpath;
addpath(genpath('src'));

config;   % load user parameters

%% 1. Load input data
BF = load_input(input);

%% 2. Build Malik-Chebyshev mesh (curvilinear) or pass through (Cartesian)
mesh_type = detect_mesh_type(BF);

if strcmp(mesh_type, 'curvilinear')

    ctype_params.structured = input.structured;

    if input.structured
        if ~isfield(ctype_params, 'n_xi') || isempty(ctype_params.n_xi)
            [n_xi_det, n_eta_det] = detect_grid(BF.X, BF.Y);
            ctype_params.n_xi  = n_xi_det;
            ctype_params.n_eta = n_eta_det;
        end
    end

    ctype_params.FD1_order = StabGrid.FD1_order;
    BF = create_ctype_mesh(BF, ctype_params);

    if ctype_params.plot
        plot_opts = struct();
        if isfield(ctype_params, 'plot_opts')
            plot_opts = ctype_params.plot_opts;
        end
        plot_ctype_BF(BF, plot_opts);
    end

else
    BF.BC_mask = false;
end

%% 3. Generate stability grid
fprintf('Generating grid for DeHNSSo ...\n');
StabGrid = generate_mesh(BF, StabGrid);
fprintf('Done!\n');

%% 5. Save output
if ~exist(output.folder, 'dir');  mkdir(output.folder);  end
filepath = fullfile(output.folder, output.filename);
save(filepath, 'StabGrid');
fprintf('Saved to %s\n', filepath);
