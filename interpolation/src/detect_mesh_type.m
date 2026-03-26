function mesh_type = detect_mesh_type(BF)
%DETECT_MESH_TYPE  Classify BF data as 'cartesian' or 'curvilinear'.
%
% For a Cartesian grid the x and y fields each contain a small set of
% unique values whose product equals (or is very close to) the total
% number of points.  A curvilinear (body-fitted) mesh has essentially
% unique (x,y) pairs at every node.
%
% Input
%   BF  - struct with fields x, y (flat 1-D column/row vectors)
%
% Output
%   mesh_type - 'cartesian' or 'curvilinear'

x = BF.X(:);
y = BF.Y(:);

% If X and Y have different lengths they are separate coordinate vectors
% (not paired points) — this is always a Cartesian grid.
if length(x) ~= length(y)
    mesh_type = 'cartesian';
    fprintf('Detected mesh type: cartesian  (%d x points, %d y points)\n', ...
        length(x), length(y));
    return
end

n     = length(x);
tol_x = max(1e-6 * (max(x) - min(x)), 1e-12);
tol_y = max(1e-6 * (max(y) - min(y)), 1e-12);

xu = uniquetol(x, tol_x, 'DataScale', 1);
yu = uniquetol(y, tol_y, 'DataScale', 1);

ratio = (length(xu) * length(yu)) / n;

if ratio > 0.9 && ratio < 1.1
    mesh_type = 'cartesian';
else
    mesh_type = 'curvilinear';
end

fprintf('Detected mesh type: %s  (%d unique x, %d unique y, %d total)\n', ...
    mesh_type, length(xu), length(yu), n);
end
