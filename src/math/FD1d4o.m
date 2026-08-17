function [ DD ] = FD1d4o( D,d )
%FD1D4O  Fourth-order first derivative on a uniform grid.
%
%   DD = FD1d4o(D, d)
%
%   Differentiates the field D along its rows, assuming uniform spacing d and
%   that the coordinate ascends as the column index decreases. Central
%   differences are used in the interior and biased five-point stencils at the
%   two boundary columns on each side. At least six columns are required.
%
%   Inputs:
%     D  — field sampled on a uniform grid (derivative taken along rows)
%     d  — grid spacing
%   Outputs:
%     DD — first derivative, same size as D

%% Check enough points for the stencil
if length(D(1,:)) < 6
    error('--- Function length too small; less than 6')
end

%% Differentiate
DD = 0*D;

% Central differences on the interior columns
DD(:,  (1+2):(end-2)) =  (1/12 * D(:,  (1+4):(end-0))...
                      -    2/3 * D(:,  (1+3):(end-1))...
                      +      0 * D(:,  (1+2):(end-2))...
                      +    2/3 * D(:,  (1+1):(end-3))...
                      -   1/12 * D(:,  (1+0):(end-4)))/ d;

% Forward differences on the two smallest-coordinate columns
DD(:,(end-1):(end-0)) =(-25/12 * D(:,(end-1):(end-0))...
                      +      4 * D(:,(end-2):(end-1))...
                      -      3 * D(:,(end-3):(end-2))...
                      +    4/3 * D(:,(end-4):(end-3))...
                      -    1/4 * D(:,(end-5):(end-4)))/ d;

% Backward differences on the two largest-coordinate columns
DD(:,  (1+0):  (1+1)) =(-25/12 * D(:,  (1+0):  (1+1))...
                      +      4 * D(:,  (1+1):  (1+2))...
                      -      3 * D(:,  (1+2):  (1+3))...
                      +    4/3 * D(:,  (1+3):  (1+4))...
                      -    1/4 * D(:,  (1+4):  (1+5)))/(-d);

end
