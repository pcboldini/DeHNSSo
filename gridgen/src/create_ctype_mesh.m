function [BF_new] = create_ctype_mesh(BF, params)
%CREATE_CTYPE_MESH  Build HNS-ready structured mesh from a C-type CSV dataset.
%
% Finds the wall geometry from the no-slip condition, builds a clean
% body-fitted target mesh (Malik-Chebyshev wall-normal, uniform streamwise),
% and interpolates the input data onto it via scatteredInterpolant.
%
% Inputs
%   BF      - struct with flat fields X, Y, U, V, W  (from load_input/CSV)
%   params  - configuration struct with fields:
%               .n_eta       : wall-normal dimension (expected wall point count)
%               .n_xi        : streamwise dimension  (expected wall point count)
%               .n_eta_new   : target wall-normal resolution (default: n_eta)
%               .n_xi_new    : target streamwise resolution  (default: n_xi)
%               .H           : (optional) wall-normal domain height
%                              (default: estimated from data)
%               .y_i         : (optional) Malik median collocation height
%                              (default: H/3)
%               .xi_range    : (optional) [frac_start frac_end] fraction of
%                              arc-length to retain (default: full range)
%               .plot        : (optional) true to plot mesh
%
% Output
%   BF_new  - struct with fields:
%               X, Y         : physical coordinates [n_eta_new x n_xi_new]
%               U, V, W      : velocity components  [n_eta_new x n_xi_new]
%               xi1D         : streamwise arc-length [1 x n_xi_new]
%               eta1D        : wall-normal Malik-Chebyshev arc-length [n_eta_new x 1]
%               nx, ny       : grid dimensions
%               H, y_i       : Malik mapping parameters
%               mesh_type    : 'curvilinear'

%% Parse parameters
structured = isfield(params, 'structured') && params.structured;

if structured
    n_eta = params.n_eta;
    n_xi  = params.n_xi;

    n_eta_new = n_eta;
    n_xi_new  = n_xi;

    %% Validate dimensions
    n_total = numel(BF.X);
    if n_eta * n_xi ~= n_total
        error('create_ctype_mesh: n_eta (%d) * n_xi (%d) = %d but data has %d points.', ...
              n_eta, n_xi, n_eta * n_xi, n_total);
    end
else
    n_eta_new = 101;
    n_xi_new  = [];
end

if isfield(params, 'n_eta_new') && ~isempty(params.n_eta_new)
    n_eta_new = params.n_eta_new;
end
if isfield(params, 'n_xi_new') && ~isempty(params.n_xi_new)
    n_xi_new = params.n_xi_new;
end

X_flat = BF.X(:);
Y_flat = BF.Y(:);
U_flat = BF.U(:);
V_flat = BF.V(:);
W_flat = BF.W(:);

%% Non-dimensionalise if reference scales are provided
if isfield(params, 'lref') && ~isempty(params.lref) && ...
   isfield(params, 'Uref') && ~isempty(params.Uref)
    fprintf('Non-dimensionalising: Uref = %.4f, lref = %.6f\n', params.Uref, params.lref);
    X_flat = X_flat / params.lref;
    Y_flat = Y_flat / params.lref;
    U_flat = U_flat / params.Uref;
    V_flat = V_flat / params.Uref;
    W_flat = W_flat / params.Uref;
end

%% Input data summary
ux = unique(X_flat);
idx_inlet = abs(X_flat - ux(1)) < 1e-10;
fprintf('Input data: %d points, %d unique x-stations\n', numel(X_flat), numel(ux));
fprintf('Inlet station (x=%.4f): %d wall-normal points, y=[%.4f, %.4f]\n', ...
        ux(1), sum(idx_inlet), min(Y_flat(idx_inlet)), max(Y_flat(idx_inlet)));
y_top = max(Y_flat(idx_inlet));
idx_top = abs(Y_flat - y_top) < 1e-10 * y_top;
ux_top = unique(X_flat(idx_top));
fprintf('Top station (y=%.4f): %d streamwise points, x=[%.4f, %.4f]\n', ...
        y_top, numel(ux_top), min(ux_top), max(ux_top));
fprintf('max U = %.4f, freestream U (y=%.4f) = %.4f\n', ...
        max(abs(U_flat)), y_top, mean(U_flat(idx_top)));

%% Step 1: Find wall points from no-slip condition (U=0, V=0)

if structured
    % Structured grid: take n_xi points with smallest |U|+|V|,
    % with minimum spatial separation to enforce one point per xi station.
    uv_flat  = abs(U_flat) + abs(V_flat);
    [~, sidx] = sort(uv_flat);
    sidx      = sidx(1 : min(2*n_xi, numel(sidx)));
    min_dist  = xi_approx_step(X_flat, Y_flat, n_xi) * 0.5;

    wall_idx = zeros(n_xi, 1);
    n_found  = 0;
    X_sel    = zeros(n_xi, 1);
    Y_sel    = zeros(n_xi, 1);

    for i = 1:length(sidx)
        if n_found == n_xi;  break;  end
        c = sidx(i);
        if n_found == 0 || min(sqrt((X_sel(1:n_found) - X_flat(c)).^2 + ...
                                    (Y_sel(1:n_found) - Y_flat(c)).^2)) > min_dist
            n_found           = n_found + 1;
            wall_idx(n_found) = c;
            X_sel(n_found)    = X_flat(c);
            Y_sel(n_found)    = Y_flat(c);
        end
    end

    wall_idx = wall_idx(1:n_found);
    n_xi     = n_found;

    if n_found ~= params.n_xi
        warning('Wall detection: found %d points but expected %d.', n_found, params.n_xi);
    end
    fprintf('Wall detection: found %d points (expected %d)\n', n_found, params.n_xi);
else
    % Unstructured scatter: wall points are where |U| + |V| ≈ 0
    uv_flat  = abs(U_flat) + abs(V_flat);
    wall_tol = 1e-6;
    wall_idx = find(uv_flat < wall_tol);
    n_xi     = numel(wall_idx);
    if isempty(n_xi_new)
        n_xi_new = n_xi;
    end
    fprintf('Wall detection (unstructured): found %d points with |U|+|V| < %g\n', n_xi, wall_tol);
end

X_wall = X_flat(wall_idx);
Y_wall = Y_flat(wall_idx);

% Plot wall detection and wait for user confirmation
vel_names  = {'U (streamwise)', 'V (wall-normal)', 'W (spanwise)'};
vel_data   = {U_flat, V_flat, W_flat};

figure;
for ip = 1:3
    subplot(4,1,ip);
    scatter(X_flat, Y_flat, 4, vel_data{ip}, 'filled');
    clim([min(vel_data{ip}) max(vel_data{ip})]);
    colormap(jet); colorbar;
    daspect([1 1 1]);
    title(vel_names{ip});
end

subplot(4,1,4);
plot(X_flat, Y_flat, 'k.', 'MarkerSize', 1);
hold on;
plot(X_wall, Y_wall, 'r.', 'MarkerSize', 6);
daspect([1 1 1]);
legend('input points', 'detected wall', 'Location', 'best');
title(sprintf('Wall detection — %d points', n_xi));
if usejava('desktop')
    input('Wall detection OK? Press Enter to continue or Ctrl+C to abort: ', 's');
end

%% Step 2: Sort wall points along the wall curve
%
% Directional nearest-neighbour traversal starting from the leftmost point.
% At each step, candidates in the backward half-plane are excluded to
% prevent zigzag near the trailing edge.

[~, idx_minx] = min(X_wall);
[~, idx_maxx] = max(X_wall);
fprintf('Wall point with min x: (x=%.4f, y=%.4f)\n', X_wall(idx_minx), Y_wall(idx_minx));
fprintf('Wall point with max x: (x=%.4f, y=%.4f)\n', X_wall(idx_maxx), Y_wall(idx_maxx));

order    = zeros(n_xi, 1);
visited  = false(n_xi, 1);
order(1) = idx_minx;
visited(idx_minx) = true;

for k = 2:n_xi
    prev = order(k-1);
    dx   = X_wall - X_wall(prev);
    dy   = Y_wall - Y_wall(prev);
    d    = dx.^2 + dy.^2;
    d(visited) = inf;

    if k > 2
        prev2   = order(k-2);
        dir_x   = X_wall(prev) - X_wall(prev2);
        dir_y   = Y_wall(prev) - Y_wall(prev2);
        d_fwd   = d;
        d_fwd(dx * dir_x + dy * dir_y < 0) = inf;
        if any(~isinf(d_fwd))
            [~, next] = min(d_fwd);
        else
            [~, next] = min(d);   % fallback: no forward candidate
        end
    else
        [~, next] = min(d);
    end

    order(k)      = next;
    visited(next) = true;
end

wall_idx = wall_idx(order);
X_wall   = X_wall(order);
Y_wall   = Y_wall(order);

% Remove coincident wall points (zero arc-length step)
dxi_wall    = sqrt(diff(X_wall).^2 + diff(Y_wall).^2);
unique_mask = [true, (dxi_wall > 0)'];
if any(~unique_mask)
    fprintf('Removing %d coincident wall points.\n', sum(~unique_mask));
    wall_idx = wall_idx(unique_mask);
    X_wall   = X_wall(unique_mask);
    Y_wall   = Y_wall(unique_mask);
    dxi_wall = sqrt(diff(X_wall).^2 + diff(Y_wall).^2);
    n_xi     = sum(unique_mask);
end

xi_wall = [0, cumsum(dxi_wall)'];   % [1 x n_xi]
fprintf('Wall arc-length: %.4f\n', xi_wall(end));

%% Step 3: Wall tangent and outward normal at each wall point
%
% Fit splines X(xi) and Y(xi) and evaluate derivatives analytically.

ppX = spline(xi_wall, X_wall);
ppY = spline(xi_wall, Y_wall);

tx = ppval(fnder(ppX), xi_wall(:))';   % [1 x n_xi]
ty = ppval(fnder(ppY), xi_wall(:))';

t_len = sqrt(tx.^2 + ty.^2);
tx = tx ./ t_len;
ty = ty ./ t_len;

nx_w = -ty;   % outward wall normal [1 x n_xi]
ny_w =  tx;

% Warn if a normal fold is detected (trailing-edge fold)
d_angle  = abs(diff(unwrap(atan2(ny_w, nx_w))));
idx_fold = find(d_angle > pi/4, 1, 'first');
if ~isempty(idx_fold)
    frac = xi_wall(idx_fold) / xi_wall(end);
    warning('Normal fold of %.1f deg at %.0f%% of arc-length. Set xi_range = [0, %.2f].', ...
        rad2deg(d_angle(idx_fold)), frac*100, frac - 0.02);
end

%% Step 4: Wall-normal domain height H

if isfield(params, 'H') && ~isempty(params.H)
    H = params.H;
    fprintf('Using prescribed H = %.5f\n', H);
else
    H_col = zeros(n_xi, 1);
    for j = 1:n_xi
        proj     = (X_flat - X_wall(j)) * nx_w(j) + (Y_flat - Y_wall(j)) * ny_w(j);
        H_col(j) = max(proj);
    end
    H = mean(H_col);
    fprintf('Estimated H = %.5f\n', H);
end

%% Step 4b: Estimate delta_99 from raw input data

eta_probe = linspace(0, H, 500)';

% Build a temporary scattered interpolant for U only
F_U_tmp = scatteredInterpolant(X_flat, Y_flat, U_flat, 'natural', 'nearest');

% Estimate delta_99 at the inflow station (first wall point)
U_ray_in = abs(F_U_tmp(X_wall(1) + eta_probe * nx_w(1), ...
                        Y_wall(1) + eta_probe * ny_w(1)));
U_edge   = U_ray_in(end);
idx99    = find(U_ray_in >= 0.99 * U_edge, 1, 'first');
if ~isempty(idx99)
    delta_99 = eta_probe(idx99);
else
    delta_99 = H;
end
fprintf('Estimated delta_99 = %.5f\n', delta_99);
fprintf('  Suggestion: H   ~ %.5f  (e.g. 4 * delta_99)\n', 4 * delta_99);
fprintf('  Suggestion: y_i ~ %.5f  (e.g. delta_99)\n',     delta_99);

if isfield(params, 'y_i') && ~isempty(params.y_i)
    y_i = params.y_i;
    fprintf('Using prescribed y_i = %.5f\n', y_i);
else
    y_i = delta_99;
    fprintf('Using y_i = delta_99 = %.5f\n', y_i);
end

%% Step 5: Malik-Chebyshev target eta coordinates

a = H * y_i / (H - 2.0 * y_i);
b = 1.0 + 2.0 * a / H;

% Descending: z(1)=+1 -> eta=H (freestream), z(end)=-1 -> eta=0 (wall)
% Matches DeHNSSo convention: freestream at index 1, wall at index end.
k_vec    = (0 : n_eta_new - 1)';
z_cheb   = cos(pi * k_vec / (n_eta_new - 1));
eta_cheb = a * (1 + z_cheb) ./ (b - z_cheb);   % [n_eta_new x 1], descending H -> 0

%% Step 6: Apply optional xi trimming

if isfield(params, 'xi_range') && ~isempty(params.xi_range)
    xi_lo = params.xi_range(1) * xi_wall(end);
    xi_hi = params.xi_range(2) * xi_wall(end);
    keep  = xi_wall >= xi_lo & xi_wall <= xi_hi;
    if sum(keep) < 2
        error('create_ctype_mesh: xi_range [%.2f, %.2f] retains fewer than 2 points.', ...
              params.xi_range(1), params.xi_range(2));
    end
    xi_wall = xi_wall(keep);
    X_wall  = X_wall(keep);
    Y_wall  = Y_wall(keep);
    nx_w    = nx_w(keep);
    ny_w    = ny_w(keep);
    fprintf('xi trimmed to %.0f%%–%.0f%% of arc-length: %d stations retained\n', ...
        params.xi_range(1)*100, params.xi_range(2)*100, sum(keep));
end

%% Step 7: Build target mesh
%
% Resample wall geometry to uniform xi, then shoot Malik-Chebyshev rays:
%   X_target(k,j) = X_wall_u(j) + eta_cheb(k) * nx_u(j)

xi1D_new = linspace(xi_wall(1), xi_wall(end), n_xi_new);   % [1 x n_xi_new]

X_wall_u = interp1(xi_wall, X_wall, xi1D_new, 'pchip');
Y_wall_u = interp1(xi_wall, Y_wall, xi1D_new, 'pchip');
nx_u     = interp1(xi_wall, nx_w,   xi1D_new, 'pchip');
ny_u     = interp1(xi_wall, ny_w,   xi1D_new, 'pchip');

% Renormalise normals after interpolation
n_len = sqrt(nx_u.^2 + ny_u.^2);
nx_u  = nx_u ./ n_len;
ny_u  = ny_u ./ n_len;

% eta_cheb [n_eta_new x 1] x [1 x n_xi_new] -> [n_eta_new x n_xi_new]
X_target = X_wall_u + eta_cheb .* nx_u;
Y_target = Y_wall_u + eta_cheb .* ny_u;

%% Step 7b: Elliptic mesh smoothing (optional)
if isfield(params, 'elliptic') && params.elliptic
    ell_opts = struct('max_iter', 5000, 'tol', 1e-6, 'omega', 1.7, 'n_uniform', 201);
    if isfield(params, 'elliptic_opts')
        flds = fieldnames(params.elliptic_opts);
        for f = 1:numel(flds), ell_opts.(flds{f}) = params.elliptic_opts.(flds{f}); end
    end
    [X_target, Y_target] = solve_elliptic_mesh(X_target, Y_target, ...
        X_wall_u, Y_wall_u, nx_u, ny_u, eta_cheb, H, ell_opts);
end

%% Step 8: Scattered interpolation of U, V, W onto target mesh

F_U = scatteredInterpolant(X_flat, Y_flat, U_flat, 'natural', 'nearest');
F_V = scatteredInterpolant(X_flat, Y_flat, V_flat, 'natural', 'nearest');
F_W = scatteredInterpolant(X_flat, Y_flat, W_flat, 'natural', 'nearest');

U_new = F_U(X_target, Y_target);
V_new = F_V(X_target, Y_target);
W_new = F_W(X_target, Y_target);

% Enforce no-slip at wall (wall is at last index: descending convention)
U_new(end,:) = 0;
V_new(end,:) = 0;
W_new(end,:) = 0;

fprintf('max U input:        %.6f\n', max(U_flat));
fprintf('max U interpolated: %.6f\n', max(U_new(:)));


%% Plot: velocities along a constant wall-normal ray

if isfield(params, 'plot_eta') && ~isempty(params.plot_eta)
    eta_ray = params.plot_eta;
else
    eta_ray = delta_99 * 5;
end
X_ray = X_wall_u + eta_ray * nx_u;
Y_ray = Y_wall_u + eta_ray * ny_u;

U_ray = F_U(X_ray(:), Y_ray(:));
V_ray = F_V(X_ray(:), Y_ray(:));
W_ray = F_W(X_ray(:), Y_ray(:));

% Raw input points near the ray: select points within a geometric band
% around the ray curve  X_ray(xi), Y_ray(xi).
ray_band = max(delta_99, eta_ray) * 0.1;
dist_sq  = inf(numel(X_flat), 1);
xi_near  = zeros(numel(X_flat), 1);
for j = 1:n_xi_new
    d2 = (X_flat - X_ray(j)).^2 + (Y_flat - Y_ray(j)).^2;
    closer = d2 < dist_sq;
    dist_sq(closer) = d2(closer);
    xi_near(closer) = xi1D_new(j);
end
near_ray = sqrt(dist_sq) < ray_band;

% FD derivatives along the ray (same scheme used downstream)
dU_ray = FD1_uneven(xi1D_new(:), U_ray(:), params.FD1_order)';
dV_ray = FD1_uneven(xi1D_new(:), V_ray(:), params.FD1_order)';
dW_ray = FD1_uneven(xi1D_new(:), W_ray(:), params.FD1_order)';

figure;
vel_ray  = {U_ray, V_ray, W_ray};
dvel_ray = {dU_ray, dV_ray, dW_ray};
raw_vel  = {U_flat(near_ray), V_flat(near_ray), W_flat(near_ray)};
raw_xi   = xi_near(near_ray);
colors   = {'b', 'r', 'k'};
labels   = {'U', 'V', 'W'};
for ip = 1:3
    subplot(3,2,2*ip-1);
    plot(raw_xi, raw_vel{ip}, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 4); hold on;
    plot(xi1D_new, vel_ray{ip}, [colors{ip} '-'], 'LineWidth', 1.2);
    ylabel(labels{ip});
    legend('raw input', 'interpolated', 'Location', 'best');
    if ip == 1; title(sprintf('Velocities at eta = %.4f (5 \\delta_{99})', eta_ray)); end

    subplot(3,2,2*ip);
    plot(xi1D_new, dvel_ray{ip}, [colors{ip} '-'], 'LineWidth', 1.2);
    ylabel(['d' labels{ip} '/d\xi']);
    if ip == 1; title('FD streamwise derivatives'); end
    if ip == 3; xlabel('\xi (arc-length)'); end
end

%% Plot: wall-normal profiles at selected xi station

if isfield(params, 'plot_xi') && ~isempty(params.plot_xi)
    pxi = params.plot_xi;
    if pxi <= 1
        j_mid = max(1, round(pxi * n_xi_new));
    else
        j_mid = round(pxi);
    end
else
    j_mid = round(n_xi_new / 2);
end
X_prof = X_target(:, j_mid);
Y_prof = Y_target(:, j_mid);

U_prof = U_new(:, j_mid);
V_prof = V_new(:, j_mid);
W_prof = W_new(:, j_mid);

% Wall-normal derivatives for diagnostic plot
use_fd_plot = isfield(params, 'plot_FD') && params.plot_FD;
if use_fd_plot
    dU_prof = FD1_uneven(eta_cheb, U_prof, params.FD1_order);
    dV_prof = FD1_uneven(eta_cheb, V_prof, params.FD1_order);
    dW_prof = FD1_uneven(eta_cheb, W_prof, params.FD1_order);
    deriv_label = 'FD wall-normal derivatives';
else
    tmp_grid.ny = n_eta_new;  tmp_grid.H = H;  tmp_grid.y_i = y_i;  tmp_grid.y = eta_cheb;
    tmp_grid = compute_cheb_operators(tmp_grid);
    dU_prof = tmp_grid.D1 * U_prof;
    dV_prof = tmp_grid.D1 * V_prof;
    dW_prof = tmp_grid.D1 * W_prof;
    deriv_label = 'Chebyshev wall-normal derivatives';
end

% Raw input points near this profile
prof_band = (xi1D_new(end) - xi1D_new(1)) / n_xi_new * 0.6;
dist_prof = sqrt((X_flat - X_wall_u(j_mid)).^2 + (Y_flat - Y_wall_u(j_mid)).^2);
eta_raw   = (X_flat - X_wall_u(j_mid)) * nx_u(j_mid) + (Y_flat - Y_wall_u(j_mid)) * ny_u(j_mid);
tau_raw   = abs((X_flat - X_wall_u(j_mid)) * (-ny_u(j_mid)) + (Y_flat - Y_wall_u(j_mid)) * nx_u(j_mid));
near_prof = tau_raw < prof_band & eta_raw >= 0 & eta_raw <= H;

figure;
vel_prof  = {U_prof, V_prof, W_prof};
dvel_prof = {dU_prof, dV_prof, dW_prof};
raw_prof  = {U_flat(near_prof), V_flat(near_prof), W_flat(near_prof)};
eta_near  = eta_raw(near_prof);
for ip = 1:3
    subplot(3,2,2*ip-1);
    plot(eta_near, raw_prof{ip}, '.', 'Color', [0.7 0.7 0.7], 'MarkerSize', 4); hold on;
    plot(eta_cheb, vel_prof{ip}, [colors{ip} '-'], 'LineWidth', 1.2);
    xlabel('\eta'); ylabel(labels{ip});
    legend('raw input', 'interpolated', 'Location', 'best');
    if ip == 1; title(sprintf('Wall-normal profiles at \\xi = %.4f (station %d/%d)', xi1D_new(j_mid), j_mid, n_xi_new)); end

    subplot(3,2,2*ip);
    plot(eta_cheb, dvel_prof{ip}, [colors{ip} '-'], 'LineWidth', 1.2);
    xlabel('\eta'); ylabel(['d' labels{ip} '/d\eta']);
    if ip == 1; title(deriv_label); end
end

%% Plot

if isfield(params, 'plot') && params.plot
    vel_names = {'U', 'V', 'W'};
    vel_new   = {U_new, V_new, W_new};

    figure;
    for ip = 1:3
        subplot(4,1,ip);
        scatter(X_target(:), Y_target(:), 4, vel_new{ip}(:), 'filled');
        clim([min(vel_new{ip}(:)) max(vel_new{ip}(:))]);
        colormap(jet); colorbar;
        daspect([1 1 1]);
        title(sprintf('Interpolated %s on target mesh', vel_names{ip}));
    end

    subplot(4,1,4);
    plot(X_flat, Y_flat, 'k.', 'MarkerSize', 1);
    hold on;
    plot(X_wall, Y_wall, 'r.', 'MarkerSize', 4);
    plot(X_target,  Y_target,  'b-', 'LineWidth', 0.5);
    plot(X_target', Y_target', 'b-', 'LineWidth', 0.5);
    daspect([1 1 1]);
    legend('input data', 'wall', 'target mesh', 'Location', 'best');
    title(sprintf('Target mesh: %d xi x %d eta', n_xi_new, n_eta_new));
end

%% Step 9: Assemble output structure

BF_new.X        = X_target;
BF_new.Y        = Y_target;
BF_new.U        = U_new;
BF_new.V        = V_new;
BF_new.W        = W_new;
BF_new.xi1D     = xi1D_new;
BF_new.eta1D    = eta_cheb;
BF_new.nx       = n_xi_new;
BF_new.ny       = n_eta_new;
BF_new.H        = H;
BF_new.y_i      = y_i;
BF_new.mesh_type = 'curvilinear';
if isfield(params, 'Re') && ~isempty(params.Re)
    BF_new.Re = params.Re;
end

fprintf('C-type mesh ready: %d streamwise x %d wall-normal  (H=%.4f, y_i=%.4f)\n', ...
    n_xi_new, n_eta_new, H, y_i);
end

function d = xi_approx_step(X_flat, Y_flat, n_xi)
Lx = range(X_flat);
Ly = range(Y_flat);
d  = sqrt(Lx^2 + Ly^2) / n_xi;
end
