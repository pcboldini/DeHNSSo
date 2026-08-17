function [ D ] = FD1d2o_uneven(x,y)
%FD1D2O_UNEVEN  First derivative of tabulated data on a non-uniform grid.
%
%   D = FD1d2o_uneven(x, y)
%
%   Computes dy/dx with a three-point finite-difference (FD) stencil that
%   accounts for the uneven spacing between grid points: central differences
%   in the interior, forward at the first node, backward at the last. Data is
%   given in columns; for a matrix, each column is differentiated on its own.
%
%   Inputs:
%     x  — grid coordinates (column vector)
%     y  — sampled data (column vector or matrix)
%   Outputs:
%     D  — dy/dx, same size as y

dx = diff(x);
D = y;

%% Central differences (interior nodes)
D(2:end-1,:) = (   dx(1:end-1,:).^2                  .* y(3:end  ,:) ...
                + (dx(2:end  ,:).^2-dx(1:end-1,:).^2).* y(2:end-1,:) ...
                -  dx(2:end  ,:).^2                  .* y(1:end-2,:)  ) ...
               ./ (dx(1:end-1,:).*dx(2:end,:).*(dx(1:end-1,:) + dx(2:end,:)));

%% Forward differences (first node)
D(1,:)       = (-  dx(1,:).^2                        .* y(3,:) ...
                + (dx(1,:)    +    dx(2,:)).^2       .* y(2,:) ...
                - (dx(2,:).^2 + 2.*dx(1,:).*dx(2,:)) .* y(1,:)  ) ...
               ./ (dx(1,:).*dx(2,:).*(dx(1,:) + dx(2,:)));

%% Backward differences (last node)
D(end,:)     = (  (dx(end-1,:).^2 + 2.*dx(end-1,:).*dx(end,:)) .* y(end,:) ...
                - (dx(end-1,:)    +    dx(end  ,:)).^2         .* y(end-1,:) ...
                +  dx(end  ,:).^2                              .* y(end-2,:)   ) ...
               ./ (dx(end-1,:).*dx(end,:).*(dx(end-1,:)+dx(end,:)));

end
