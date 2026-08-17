function [ynd, D1y, D2y] = MappingMalik(y_max, y_i, eta, D1eta, D2eta)
%MAPPINGMALIK  Map a Chebyshev grid [-1,1] to a clustered physical grid [0,y_max].
%
%   [ynd, D1y, D2y] = MappingMalik(y_max, y_i, eta, D1eta, D2eta)
%
%   Malik's algebraic mapping, with half the nodes below the median height y_i.
%   Also returns the physical differentiation matrices for the mapped grid.
%
%   Inputs:
%     y_max         — domain height
%     y_i           — median height (half the nodes lie in [0, y_i])
%     eta           — Chebyshev nodes on [-1, 1]
%     D1eta, D2eta  — Chebyshev differentiation matrices
%   Outputs:
%     ynd           — physical wall-normal grid [0, y_max]
%     D1y, D2y      — first / second physical differentiation matrices

ynd = (y_i*y_max*(1+eta) ./ (y_max - eta*(y_max - 2*y_i)))';

if nargout > 1
    detady   =  (y_max - eta*(y_max - 2*y_i)).^2 / (2*y_i*y_max*(y_max - y_i));
    d2etady2 = -2*detady.^2*(y_max - 2*y_i) ./ (y_max - eta*(y_max - 2*y_i));
    D1y = diag(detady)*D1eta;
    D2y = diag(d2etady2)*D1eta + diag(detady.^2)*D2eta;
end
end
