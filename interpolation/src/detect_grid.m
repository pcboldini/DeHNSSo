function [n_streamwise, n_wallnormal] = detect_grid(x, y)
%DETECT_GRID  Detect structured grid dimensions from a flat ordered point list.
%
% Uses the physical aspect ratio of the point cloud to select the factor
% pair (n_wallnormal, n_streamwise) whose ratio best matches Lx/Ly.
% Works well for curvilinear meshes where the domain extent roughly
% reflects the resolution ratio (e.g. C-type airfoil meshes).
%
% Inputs
%   x, y          - flat 1-D coordinate vectors (ordered, not scattered)
%
% Outputs
%   n_streamwise  - slow / outer index  (number of xi stations)
%   n_wallnormal  - fast / inner index  (number of eta levels)

x = x(:);
y = y(:);
N = length(x);

Lx     = range(x);
Ly     = range(y);
aspect = Lx / Ly;

divs = find(rem(N, 1:N) == 0);

best_err     = inf;
n_streamwise = N;
n_wallnormal = 1;

for i = 1:length(divs)
    nw = divs(i);
    ns = N / nw;

    % Require both dimensions to have at least 5 points
    if nw < 5 || ns < 5
        continue;
    end

    % Streamwise direction is assumed to be the longer one
    if ns >= nw
        err = abs(ns/nw - aspect);
        if err < best_err
            best_err     = err;
            n_streamwise = ns;
            n_wallnormal = nw;
        end
    end
end

fprintf('Detected grid: %d streamwise x %d wall-normal\n', n_streamwise, n_wallnormal);
end
