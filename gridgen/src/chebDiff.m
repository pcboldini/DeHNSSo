function D = chebDiff(N)
%CHEBDIFF  Chebyshev spectral differentiation matrix.
%
%   D = chebDiff(N) returns the N x N first-order differentiation matrix
%   for Chebyshev points x_k = -cos(pi*k/(N-1)), k = 0..N-1.
%   Points are ordered from -1 to +1.
%
%   Reference: Trefethen, Spectral Methods in MATLAB, SIAM (2000), cheb.m

if N == 0; D = 0; return; end

N = N - 1;
x = -cos(pi*(0:N)/N)';                          % Chebyshev points [-1, +1]
c = [2; ones(N-1,1); 2] .* (-1).^(0:N)';        % barycentric weights
X  = repmat(x, 1, N+1);
dX = X - X';                                     % pairwise differences
D  = (c * (1./c)') ./ (dX + eye(N+1));           % off-diagonal entries
D  = D - diag(sum(D, 2));                         % diagonal: negative row sum

end
