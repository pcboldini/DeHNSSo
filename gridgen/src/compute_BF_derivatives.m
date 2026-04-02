function Grid = compute_BF_derivatives(Grid)
%COMPUTE_BF_DERIVATIVES  First derivatives of base-flow velocities U, V, W.
%
%   Wall-normal (eta):  dyF(:,i) = D1 * F(:,i)   [spectral, Malik-Chebyshev]
%   Streamwise  (xi):   dxF(j,:) via FD on xi1D  [finite difference]

Grid.dyU = Grid.D1 * Grid.U;
Grid.dyV = Grid.D1 * Grid.V;
Grid.dyW = Grid.D1 * Grid.W;

xi = Grid.xi1D(:);
Grid.dxU = zeros(Grid.ny, Grid.nx);
Grid.dxV = zeros(Grid.ny, Grid.nx);
Grid.dxW = zeros(Grid.ny, Grid.nx);
for j = 1:Grid.ny
    Grid.dxU(j,:) = FD1_uneven(xi, Grid.U(j,:)', Grid.FD1_order)';
    Grid.dxV(j,:) = FD1_uneven(xi, Grid.V(j,:)', Grid.FD1_order)';
    Grid.dxW(j,:) = FD1_uneven(xi, Grid.W(j,:)', Grid.FD1_order)';
end
end
