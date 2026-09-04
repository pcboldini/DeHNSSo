function StabGrid = main_gridgen(input, params, output)
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
% Build a DeHNSSo-ready StabGrid from a raw base-flow file.
%
%   StabGrid = main_gridgen(input, params, output)
%
%   Three input paths are supported, auto-selected by classify_input:
%     'cartesian'              — 1-D X, 1-D Y axes (e.g. bf_blasius.mat)
%     'curvilinear_structured' — 2-D X, 2-D Y body-fitted grid
%                                (e.g. bf_sweptwing_flat.mat, bf_sweptwing_hump.mat)
%     'unstructured'           — flat scatter cloud (e.g. m3j.csv)
%
%   input   struct with fields:
%     .folder     directory containing the BF file
%     .filename   BF filename
%     .format     'mat' | 'bin' | 'csv'
%     .structured (optional) — hint passed to classify_input
%
%   params  struct with the union of fields consumed by the chosen resampler
%           (n_eta_new / n_xi_new / y_i / H / xi_range / xi_trim_* / Re /
%            lref / Uref / nu / rescale / elliptic / plot / plot_opts ...)
%           plus grid-generation knobs (FD_xi_order_1, FD_xi_order_2, FD_eta_method,
%           FD_eta_order). See the per-case callers in gridgen/benchmark/.
%
%   output  struct with fields:
%     .folder     directory where the output .mat is written
%     .filename   output filename (e.g. 'DeHNSSo_input_hump.mat')
%
%   Returns the StabGrid struct and writes it to disk.

%% Inputs
% Name   size         unit   explanation
% input  (1,1) struct [-]    base-flow source: folder, filename, format ('mat', 'bin' or 'csv') and an optional structured hint
% params (1,1) struct [-]    grid-generation parameters
% output (1,1) struct [-]    where to write the generated grid: folder and filename

%% Outputs
% Name     size         unit   explanation
% StabGrid (1,1) struct [-]    stability grid, coordinate metrics and base-flow fields on it

%% 0. Banner
fprintf('\n');
fprintf(' ______     _   _  _   _  _____ _____       \n');
fprintf(' |  _  \\   | | | || \\ | |/  ___/  ___|      \n');
fprintf(' | | | |___| |_| ||  \\| |\\ `--.\\ `--.  ___  \n');
fprintf(' | | | / _ \\  _  || . ` | `--. \\`--. \\/ _ \\ \n');
fprintf(' | |/ /  __/ | | || |\\  |/\\__/ /\\__/ / (_) |\n');
fprintf(' |___/ \\___\\_| |_/\\_| \\_/\\____/\\____/ \\___/ \n');
fprintf('                                       v1.2  \n');
fprintf('  Delft Harmonic Navier-Stokes Solver\n');
fprintf('  =========================================\n');

%% 1. Load raw BF
BF = load_bf(input);

%% 1b. Wall forcing: size the clustering to params.forcing (see prepare_forcing)
% params.forcing is one spec, a cell of specs, or an actuator with vector xc.
forcing = {};
if isfield(params, 'forcing') && ~isempty(params.forcing)
    specs        = flatten_forcing(params.forcing);
    forcing      = cellfun(@prepare_forcing, specs, 'UniformOutput', false);
    params.mug   = cellfun(@(f) f.xc, forcing);
    if ~isfield(params, 'spread') || isempty(params.spread)   % default: refine to the forcing width
        params.spread = cellfun(@(f) f.w, forcing);
    end
    fprintf('Wall forcing: %d source(s) -> refining mesh\n', numel(forcing));
end

% Clustering is set through mug (centre), spread (width) and refine (strength).
if isfield(params, 'sig')
    error('main_gridgen:legacySig', ...
          ['params.sig is no longer an input: use params.spread, the half-width ' ...
           'the refinement reaches, in physical x (sig = 3*spread/half-domain).']);
end
if isfield(params, 'ag')
    error('main_gridgen:legacyAg', ...
          ['params.ag is no longer an input: use params.refine = 1/ag, how many ' ...
           'times finer dxi is at the centre than on an equidistant mesh.']);
end
% refine says how many times finer dxi is at the centre than on an equidistant
% mesh.
if isfield(params, 'refine') && ~isempty(params.refine)
    params.ag = 1 ./ params.refine;
end

%% 2. Classify and dispatch to the appropriate resampler
class = classify_input(BF, input);

if isfield(params,'plot') && params.plot
    plot_raw_bf(BF, class);
end

switch class
    case 'cartesian'
        % resample_cartesian handles both gridded and flat-scatter inputs,
        % and defaults to source resolution when n_eta_new / n_xi_new are
        % empty (effective no-op resample → just canonicalises shape).
        BF_rs = resample_cartesian(BF, params);

    case 'curvilinear_structured'
        BF_rs = resample_curvilinear(BF, params);

    case 'unstructured'
        BF_rs = resample_unstructured(BF, params);

    otherwise
        error('gridgen: unknown input class ''%s''.', class);
end

%% 3. Guarantee Re/lref/Uref/nu are on BF_rs before build_stab_grid sees it
BF_rs = ensure_metadata(BF, BF_rs, params);

%% 4. Optional plot of the resampled BF
if isfield(params, 'plot') && params.plot
    plot_opts = struct();
    if isfield(params, 'plot_opts');  plot_opts = params.plot_opts;  end
    plot_resampled_bf(BF_rs, plot_opts);
end

%% 5. Build the DeHNSSo StabGrid
fprintf('Generating StabGrid ...\n');
StabGrid = build_stab_grid(BF_rs, params_to_grid(params));
fprintf('StabGrid ready.\n');

%% 5b. Store the forcing profile on the grid, aligned to xun for the caller
% Sample at the PHYSICAL wall-x per station (xun itself on a flat plate, the wall
% row of the 2-D X on a curved wall), so the profile is correct on any mesh.
if ~isempty(forcing)
    xun = StabGrid.xun(1,:);
    if isfield(StabGrid, 'X') && ~isvector(StabGrid.X)
        [~, wr] = min(abs(StabGrid.y(:)));      % wall row (y ~ 0)
        x_phys  = StabGrid.X(wr, :);
    else
        x_phys  = xun;                          % cartesian: xun is physical x
    end
    F = struct('xc', {}, 'w', {}, 'x', {}, 'v', {}, 'kind', {});
    for k = 1:numel(forcing)
        f    = forcing{k};
        F(k) = struct('xc', f.xc, 'w', f.w, 'x', xun, 'v', f.eval(x_phys), 'kind', f.kind);
    end
    StabGrid.forcing = F;
    fprintf('Stored %d forcing profile(s) on %d grid stations\n', numel(F), numel(xun));
end

%% 6. Write to disk
if ~exist(output.folder, 'dir');  mkdir(output.folder);  end
filepath = fullfile(output.folder, output.filename);
save(filepath, 'StabGrid');
fprintf('Saved to %s\n', filepath);

end


function Grid = params_to_grid(params)
% Extract the grid-generation subset of params into a struct that
% build_stab_grid understands.
Grid = struct();
passthrough = {'Re','plot','FD_xi_order_1','FD_xi_order_2','FD_eta_method','FD_eta_order'};
for k = 1:numel(passthrough)
    f = passthrough{k};
    if isfield(params, f) && ~isempty(params.(f))
        Grid.(f) = params.(f);
    end
end
% Sensible defaults when the caller doesn't set them
if ~isfield(Grid, 'FD_xi_order_1');   Grid.FD_xi_order_1   = 4;      end
if ~isfield(Grid, 'FD_xi_order_2');   Grid.FD_xi_order_2   = 2;      end
if ~isfield(Grid, 'FD_eta_method');  Grid.FD_eta_method  = 'cheb'; end
if ~isfield(Grid, 'FD_eta_order');Grid.FD_eta_order= 4;      end
if ~isfield(Grid, 'plot');        Grid.plot        = false;  end
end
