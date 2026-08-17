function DD = FD1d4o(D, d)
%FD1D4O  Fourth-order first derivative on a uniform grid.
%
%   DD = FD1d4o(D, d)
%
%   Differentiates D along its rows (uniform spacing d), taking the coordinate
%   to ascend as the column index decreases. Central 5-point differences in the
%   interior, fourth-order one-sided differences at the two end columns.
%
%   Inputs:
%     D  — data, differentiated along its rows
%     d  — grid spacing
%   Outputs:
%     DD — first derivative, same size as D

if size(D, 2) < 6
    error('FD1d4o: need at least 6 columns.')
end

DD = 0*D;

% Central differences (interior columns)
DD(:, 3:end-2) = ( 1/12*D(:, 5:end)   - 2/3*D(:, 4:end-1) ...
                 + 2/3*D(:, 2:end-3)  - 1/12*D(:, 1:end-4) ) / d;

% Fourth-order one-sided differences at the two end columns
DD(:, end-1:end) = (-25/12*D(:, end-1:end) + 4*D(:, end-2:end-1) ...
                   -    3*D(:, end-3:end-2) + 4/3*D(:, end-4:end-3) ...
                   -  1/4*D(:, end-5:end-4)) / d;
DD(:, 1:2)       = (-25/12*D(:, 1:2)       + 4*D(:, 2:3) ...
                   -    3*D(:, 3:4)         + 4/3*D(:, 4:5) ...
                   -  1/4*D(:, 5:6)) / (-d);
end
