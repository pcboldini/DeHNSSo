function [ ynd,D1y,D2y ] = MappingMalik( y_max,y_i,eta,D1eta,D2eta )
%MAPPINGMALIK  Map the Chebyshev grid onto the physical wall-normal domain.
%
%   [ynd, D1y, D2y] = MappingMalik(y_max, y_i, eta, D1eta, D2eta)
%
%   Maps the Chebyshev grid on [-1, 0, 1] to the physical domain
%   [0, y_i, y_max] with the algebraic (Malik) stretching, and transforms the
%   pseudo-spectral differentiation matrices onto that grid.
%
%   Inputs:
%     y_max  — outer domain height
%     y_i    — half of the nodes are clustered below this height
%     eta    — Chebyshev grid on [-1, 1]
%     D1eta  — first-derivative Chebyshev matrix
%     D2eta  — second-derivative Chebyshev matrix
%   Outputs:
%     ynd    — physical (non-dimensional) grid
%     D1y    — first-derivative matrix on the physical grid
%     D2y    — second-derivative matrix on the physical grid

% Transform Chebyshev to physical grid
ynd = (y_i*y_max*(1+eta)./(y_max - eta*(y_max - 2*y_i)))';

if nargout > 1
    % Chain-rule scale factors for the coordinate transform
    detady   = (y_max - eta*(y_max - 2*y_i)).^2/(2*y_i*y_max*(y_max - y_i));
    d2etady2 = -2*detady.^2*(y_max - 2*y_i)./(y_max - eta*(y_max - 2*y_i));

    D1y = diag(detady)*D1eta;    
    D2y = diag(d2etady2)*D1eta + diag(detady.^2)*D2eta;
end

end
