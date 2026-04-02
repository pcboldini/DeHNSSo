function plot_ctype_BF(BF, opts)
%PLOT_CTYPE_BF  Visualise the interpolated base-flow on the new curvilinear mesh.
%
% Produces two figures:
%   1. Filled contour of streamwise velocity U on the physical (X,Y) mesh
%   2. Wall-normal profiles of U at several xi stations
%
% Optional opts fields:
%   .xi_slices   : xi indices or fractions [0,1] to plot wall-normal profiles
%                  (default: 8 equally spaced stations)
%   .eta_slices  : eta indices or fractions [0,1] to plot streamwise profiles

if nargin < 2;  opts = struct();  end

X = BF.X;   % [n_eta x n_xi]
Y = BF.Y;
U = BF.U;

[n_eta, n_xi] = size(U);

%% Figure 1: U contour on physical mesh
figure;
contourf(X, Y, U, 40, 'LineColor', 'none');
colorbar;
daspect([1 1 1]);
xlabel('x');  ylabel('y');
title('Streamwise velocity U');

%% Figure 2: Wall-normal profiles at several xi stations
if isfield(opts, 'xi_slices') && ~isempty(opts.xi_slices)
    xi_slices = opts.xi_slices;
    if all(xi_slices <= 1) && all(xi_slices >= 0)
        slice_idx = max(1, round(xi_slices * n_xi));
    else
        slice_idx = xi_slices;
    end
else
    slice_idx = round(linspace(1, n_xi, 8));
end

figure;  hold on;
colors = lines(numel(slice_idx));
for k = 1:numel(slice_idx)
    j = slice_idx(k);
    plot(U(:, j), BF.eta1D, '-', 'Color', colors(k,:), ...
         'DisplayName', sprintf('\\xi = %.3f', BF.xi1D(j)));
end
xlabel('U');  ylabel('\eta');
title('Wall-normal profiles of U');
legend('Location', 'best');
grid on;

%% Figure 3: Streamwise profiles at several eta stations
if isfield(opts, 'eta_slices') && ~isempty(opts.eta_slices)
    eta_slices = opts.eta_slices;
    if all(eta_slices <= 1) && all(eta_slices >= 0)
        eta_idx = max(1, round(eta_slices * n_eta));
    else
        eta_idx = eta_slices;
    end

    figure;  hold on;
    colors = lines(numel(eta_idx));
    for k = 1:numel(eta_idx)
        i = eta_idx(k);
        plot(BF.xi1D, U(i, :), '-', 'Color', colors(k,:), ...
             'DisplayName', sprintf('\\eta = %.3f', BF.eta1D(i)));
    end
    xlabel('\xi');  ylabel('U');
    title('Streamwise profiles of U');
    legend('Location', 'best');
    grid on;
end

end
