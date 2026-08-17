function StabGrid_cheap = make_cheap_stabgrid(StabGrid)
%MAKE_CHEAP_STABGRID  Build a lower-order copy of StabGrid for the
%   preconditioner.
%
%   StabGrid_cheap = make_cheap_stabgrid(StabGrid)
%
%   Builds a lower-order copy of StabGrid, used only to assemble a cheaper
%   preconditioner Jacobian (the Liu/Brehm cheap-stencil approach).
%
%   The residual is still computed with the full-order operators (4th-order
%   finite differences, FD, in eta and xi, from the original StabGrid). The
%   Newton preconditioner instead uses operators from this cheap-stencil
%   StabGrid:
%     D1, D2  : 2nd-order FD in eta (3-point central, sparser than 4th-order)
%     D_xi    : 2nd-order central FD in xi (3-point central, sparser than
%                the 4th-order 5-point Fornberg stencil)
%
%   Because the stencils are sparser, the resulting per-mode LHS matrices
%   have roughly half the nonzero count of the full-order version, so the
%   per-mode LU factor is correspondingly lighter too.
%
%   Reference:
%     Bi & Liu 2024 (START_LNS), Algorithm 1 (NAS-FGMRES): compute the
%     residual at full accuracy, but assemble the Jacobian with cheaper
%     stencils.
%     Haas, Browne, Fasel & Brehm 2020 (JCP): approximate-Jacobian
%     preconditioner for the time-spectral linearized Navier-Stokes solver.

StabGrid_cheap = StabGrid;

%% 2nd-order FD in eta (only valid when the original was FD; assert this)
assert(isfield(StabGrid,'FD_eta_method') && strcmpi(StabGrid.FD_eta_method, 'fd'), ...
    ['make_cheap_stabgrid expects the source grid to use FD in eta ' ...
     '(StabGrid.FD_eta_method = ''fd'').  Chebyshev grids are already dense ' ...
     'in eta and cannot be made "cheaper" by lowering order.']);

% compute_fd_operators lives in the grid-generation module and isn't on the
% path when solve_newton runs from src/; add it on demand. The module is
% gridgen/src in world_git (renamed from Interpol/src in earlier trees).
if exist('compute_fd_operators','file') ~= 2
    here = fileparts(mfilename('fullpath'));
    cands = { fullfile(here, '..', '..', 'gridgen',  'src'), ...
              fullfile(here, '..', '..', 'Interpol', 'src') };
    found = false;
    for k = 1:numel(cands)
        if exist(cands{k}, 'dir'); addpath(cands{k}); found = true; break; end
    end
    if ~found
        error('make_cheap_stabgrid:noFDops', ...
            'compute_fd_operators not found in gridgen/src or Interpol/src under %s', ...
            fullfile(here, '..', '..'));
    end
end

StabGrid_cheap.FD_eta_order = 2;
StabGrid_cheap = compute_fd_operators(StabGrid_cheap);   % rebuilds D1, D2 at order 2

%% 2nd-order central FD in xi (sparse banded, bandwidth 3)
x = StabGrid.xi1D(:);
n = numel(x);
rows = zeros(3*n, 1);
cols = zeros(3*n, 1);
vals = zeros(3*n, 1);
k = 0;
for i = 1:n
    if i == 1
        % 2nd-order forward difference (-3/2, 2, -1/2 scaled for non-uniform
        % spacing), used at the leading two stations. Fine here since this
        % operator is only used for the preconditioner, not the residual.
        h1 = x(2) - x(1);
        h2 = x(3) - x(1);
        cols(k+1:k+3) = [1, 2, 3];
        rows(k+1:k+3) = [1, 1, 1];
        % Coefficients for f'(x1) using x1,x2,x3
        vals(k+1) = -(h1 + h2) / (h1 * h2);
        vals(k+2) =  h2 / (h1 * (h2 - h1));
        vals(k+3) = -h1 / (h2 * (h2 - h1));
        k = k + 3;
    elseif i == n
        h1 = x(n)   - x(n-1);
        h2 = x(n)   - x(n-2);
        cols(k+1:k+3) = [n-2, n-1, n];
        rows(k+1:k+3) = [n,   n,   n];
        % Coefficients for f'(xn) using xn-2,xn-1,xn
        vals(k+1) =  h1 / (h2 * (h2 - h1));
        vals(k+2) = -h2 / (h1 * (h2 - h1));
        vals(k+3) =  (h1 + h2) / (h1 * h2);
        k = k + 3;
    else
        % 2nd-order central on a non-uniform 3-point stencil
        h_m = x(i)   - x(i-1);
        h_p = x(i+1) - x(i);
        cols(k+1:k+3) = [i-1, i, i+1];
        rows(k+1:k+3) = [i,   i, i];
        vals(k+1) = -h_p / (h_m * (h_m + h_p));
        vals(k+2) =  (h_p - h_m) / (h_m * h_p);
        vals(k+3) =  h_m / (h_p * (h_m + h_p));
        k = k + 3;
    end
end
StabGrid_cheap.D_xi = sparse(rows(1:k), cols(1:k), vals(1:k), n, n);

%% Recompute dxi (scalar, used by some downstream paths)
StabGrid_cheap.dxi = x(2) - x(1);
end
