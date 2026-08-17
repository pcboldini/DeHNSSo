function Grid = compute_bf_derivatives(Grid)
%COMPUTE_BF_DERIVATIVES  First (and second) physical derivatives of base-flow fields.
%
%   Always differentiates the velocity components U, V, W and stores
%   ∂U/∂x, ∂U/∂y, … as Grid.dxU, Grid.dyU, etc.
%   If Grid.rho / Grid.T are present (compressible base flow), it also
%   differentiates them and stores Grid.dxrho/dyrho, Grid.dxT/dyT.
%   These are the physical-space derivatives the HNS LHS assembly expects,
%   not the raw computational-space derivatives ∂·/∂ξ, ∂·/∂η.
%
%   Second physical derivatives ∂²/∂x², ∂²/∂y², ∂²/∂x∂y of U,V,W,T (compressible)
%   are also computed (Grid.dxxU, Grid.ddyU, Grid.dxyU/dxyV, …) — required by
%   LHS_compressible for the viscous terms — via the FULL curvilinear chain rule
%   using first AND second inverse metrics (xixx, etaxx, xixy, etaxy, …):
%       ∂²·/∂x² = ξ_x²·∂_ξξ + 2ξ_xη_x·∂_ξη + η_x²·∂_ηη + ξ_xx·∂_ξ + η_xx·∂_η
%   and analogously for ∂²/∂y² and ∂²/∂x∂y.  Reduces to ∂²/∂ξ², ∂²/∂η², ∂²/∂ξ∂η
%   on Cartesian grids (off-diagonal and second metrics zero).
%
%   The first-derivative chain rule uses the inverse metric tensors set earlier
%   in build_stab_grid:
%       ∂·/∂x = ξ_x · ∂·/∂ξ + η_x · ∂·/∂η
%       ∂·/∂y = ξ_y · ∂·/∂ξ + η_y · ∂·/∂η
%   On Cartesian grids (ξ_x = η_y = 1, ξ_y = η_x = 0) this reduces to the
%   identity: ∂·/∂x = ∂·/∂ξ, ∂·/∂y = ∂·/∂η.
%
%   Discretisation:
%     ∂/∂η : D1_bf — Grid.D1 in the interior, local one-sided Fornberg
%            stencils at the wall and freestream rows. The dense Chebyshev
%            boundary rows in Grid.D1 are necessary for the stability
%            operator (spectral accuracy at the (0,0) MFD mode), but they
%            ring on tiny non-smoothness in interpolated CFD base flows
%            (linear scattered interp has piecewise-constant slope) and
%            blow up the wall shear by O(N²). A local FD stencil avoids
%            this without affecting the stability operator.
%     ∂/∂ξ : deriv1_uneven on the (possibly non-uniform) Grid.xi1D
%     ∂²/∂ξ²: deriv2_uneven (order Grid.FD_xi_order_2, default 4)

% ---- Boundary-safe D1 for base-flow differentiation -----------------------
N     = Grid.ny;
y     = Grid.y(:);
p_bd  = 4;                              % 5-point one-sided stencil at boundaries
D1_bf = Grid.D1;
% Wall row (i=N): use one-sided stencil from the wall up
idx_w = (N - p_bd) : N;
w_w   = fornberg_weights_loc(y(N), y(idx_w), 1);
D1_bf(N, :)     = 0;
D1_bf(N, idx_w) = w_w;
% Freestream row (i=1): one-sided stencil from the top down
idx_t = 1 : (1 + p_bd);
w_t   = fornberg_weights_loc(y(1), y(idx_t), 1);
D1_bf(1, :)     = 0;
D1_bf(1, idx_t) = w_t;

% ---- Differentiate each base-flow field present on Grid -------------------
fields = {'U','V','W'};
if isfield(Grid, 'rho'); fields{end+1} = 'rho'; end
if isfield(Grid, 'T');   fields{end+1} = 'T';   end

xi             = Grid.xi1D(:);
FD_xi_order_2  = getfielddef(Grid, 'FD_xi_order_2', 4);   % 4th-order matches heritage FD1d4o∘FD1d4o
% Second base-flow derivatives are only needed by LHS_compressible;
% skip them entirely for incompressible BF (no ρ, no T on Grid).
is_compr       = isfield(Grid, 'rho') && isfield(Grid, 'T');

n_provided = 0;
for k = 1:numel(fields)
    fname = fields{k};
    F     = Grid.(fname);

    has_dx  = isfield(Grid, ['dx' fname]) && ~isempty(Grid.(['dx' fname]));
    has_dy  = isfield(Grid, ['dy' fname]) && ~isempty(Grid.(['dy' fname]));
    have_d1 = has_dx && has_dy;
    if xor(has_dx, has_dy)
        given = 'dx'; if has_dy; given = 'dy'; end
        warning('compute_bf_derivatives:PartialDeriv', ...
            'Base-flow field %s: %s%s provided without its pair; both dx%s and dy%s recomputed by FD.', ...
            fname, given, fname, fname, fname);
    end
    if have_d1
        n_provided = n_provided + 2;
    else
        dF_deta = D1_bf * F;
        dF_dxi  = zeros(Grid.ny, Grid.nx);
        for j = 1:Grid.ny
            dF_dxi(j,:) = deriv1_uneven(xi, F(j,:)', Grid.FD_xi_order_1)';
        end
        Grid.(['dx' fname]) = Grid.xix .* dF_dxi + Grid.etax .* dF_deta;
        Grid.(['dy' fname]) = Grid.xiy .* dF_dxi + Grid.etay .* dF_deta;
    end
end

% Second physical derivatives needed by the compressible viscous terms, via the
% general curvilinear chain rule (reduces to the Cartesian identity when the
% off-diagonal/second metrics vanish).  Precomputed here so the solver never
% differentiates the base flow; read as Grid.dxxU/ddyU/dxyU/dxyV by
% LHS_compressible(_curvi) and CLST.
%   d_xx F = xix^2 F_xixi + 2 xix etax F_xieta + etax^2 F_etaeta + xixx F_xi + etaxx F_eta
%   d_yy F = xiy^2 F_xixi + 2 xiy etay F_xieta + etay^2 F_etaeta + xiyy F_xi + etayy F_eta
%   d_xy F = xix xiy F_xixi + (xix etay + xiy etax) F_xieta + etax etay F_etaeta
%            + xixy F_xi + etaxy F_eta
if is_compr
    for nm = {'U','V','W','T'}
        fn = nm{1};
        F  = Grid.(fn);
        % computational derivatives of the base-flow field
        F_eta    = D1_bf * F;                          % d/deta
        F_etaeta = Grid.D2 * F;                        % d2/deta2
        F_xi     = zeros(Grid.ny, Grid.nx);
        F_xixi   = zeros(Grid.ny, Grid.nx);
        for j = 1:Grid.ny
            F_xi(j,:)   = deriv1_uneven(xi, F(j,:)', Grid.FD_xi_order_1)';
            F_xixi(j,:) = deriv2_uneven(xi, F(j,:)', FD_xi_order_2)';
        end
        F_xieta = D1_bf * F_xi;                         % d2/dxi deta (= d/deta of d/dxi)
        % physical second derivatives (chain rule)
        provided = @(f) isfield(Grid,f) && ~isempty(Grid.(f));
        if provided(['dxx' fn]); n_provided=n_provided+1; else
            Grid.(['dxx' fn]) = Grid.xix.^2.*F_xixi + 2*Grid.xix.*Grid.etax.*F_xieta ...
                + Grid.etax.^2.*F_etaeta + Grid.xixx.*F_xi + Grid.etaxx.*F_eta;
        end
        if provided(['ddy' fn]); n_provided=n_provided+1; else
            Grid.(['ddy' fn]) = Grid.xiy.^2.*F_xixi + 2*Grid.xiy.*Grid.etay.*F_xieta ...
                + Grid.etay.^2.*F_etaeta + Grid.xiyy.*F_xi + Grid.etayy.*F_eta;
        end
        if any(strcmp(fn,{'U','V'}))
            if provided(['dxy' fn]); n_provided=n_provided+1; else
                Grid.(['dxy' fn]) = Grid.xix.*Grid.xiy.*F_xixi ...
                    + (Grid.xix.*Grid.etay + Grid.xiy.*Grid.etax).*F_xieta ...
                    + Grid.etax.*Grid.etay.*F_etaeta + Grid.xixy.*F_xi + Grid.etaxy.*F_eta;
            end
        end
    end
end

if n_provided > 0
    fprintf('  compute_bf_derivatives: using %d provided base-flow derivative field(s)\n', n_provided);
end
end


function v = getfielddef(s, name, default)
if isfield(s, name) && ~isempty(s.(name))
    v = s.(name);
else
    v = default;
end
end


function w = fornberg_weights_loc(z, x, m)
% Local copy of the Fornberg weights helper from compute_fd_operators.m.
% Returns a row vector of weights for the m-th derivative at z on stencil x.
n  = numel(x) - 1;
c  = zeros(n+1, m+1);
c1 = 1;
c4 = x(1) - z;
c(1,1) = 1;

for i = 1:n
    mn  = min(i, m);
    c2  = 1;
    c5  = c4;
    c4  = x(i+1) - z;

    for j = 0:i-1
        c3 = x(i+1) - x(j+1);
        c2 = c2 * c3;

        for k = mn:-1:1
            c(i+1, k+1) = c1 * (k * c(i, k) - c5 * c(i, k+1)) / c2;
        end
        c(i+1, 1) = -c1 * c5 * c(i, 1) / c2;

        for k = mn:-1:1
            c(j+1, k+1) = (c4 * c(j+1, k+1) - k * c(j+1, k)) / c3;
        end
        c(j+1, 1) = c4 * c(j+1, 1) / c3;
    end

    c1 = c2;
end

w = c(:, m+1)';
end
