function class = classify_input(BF, input)
%CLASSIFY_INPUT  Decide which resampling path a BF struct should take.
%
%   class = classify_input(BF, input)
%
%   'cartesian' is a structured rectangular grid, X and Y 1-D vectors with U/V/W
%   2-D; 'curvilinear_structured' is body-fitted, X and Y both 2-D (ny × nx) with
%   Y varying along both axes; 'unstructured' is a flat scatter of paired points
%   with no implied grid, where the wall is not known a priori.
%
%   input.structured is taken as a strong hint; otherwise the shape and
%   uniqueness of BF.X, BF.Y decide.

use_structured_hint = isfield(input, 'structured') && ~isempty(input.structured);
structured_hint     = use_structured_hint && input.structured;

X = BF.X;
Y = BF.Y;

% Vector inputs with different lengths → Cartesian (separate x, y axes).
nx_inferred = NaN; ny_inferred = NaN;
if isvector(X) && isvector(Y) && numel(X) ~= numel(Y)
    class = 'cartesian';
    nx_inferred = numel(X);
    ny_inferred = numel(Y);

% 2-D paired arrays → body-fitted structured (curvilinear) OR Cartesian
% meshgrid. Disambiguate by checking whether Y is constant along each row.
elseif ismatrix(X) && ismatrix(Y) && ~isvector(X) && ~isvector(Y) && isequal(size(X), size(Y))
    sz = size(X);   % typically [ny, nx]
    ny_inferred = sz(1);
    nx_inferred = sz(2);
    if y_is_row_constant(Y)
        class = 'cartesian';
    else
        class = 'curvilinear_structured';
    end

% Paired flat vectors. Distinguish structured (curvilinear on scatter that
% happens to be column-major) from unstructured (truly scattered points).
elseif isvector(X) && isvector(Y) && numel(X) == numel(Y)
    if structured_hint
        class = 'curvilinear_structured';
    else
        n    = numel(X);
        tolX = max(1e-6 * (max(X) - min(X)), 1e-12);
        tolY = max(1e-6 * (max(Y) - min(Y)), 1e-12);
        nu_x = numel(uniquetol(X, tolX, 'DataScale', 1));
        nu_y = numel(uniquetol(Y, tolY, 'DataScale', 1));
        ratio = (nu_x * nu_y) / n;
        if ratio > 0.9 && ratio < 1.1
            class = 'cartesian';
            nx_inferred = nu_x;
            ny_inferred = nu_y;
        else
            class = 'unstructured';
        end
    end

else
    error('classify_input:shape', ...
        'Unsupported BF.X / BF.Y shapes: X=%s, Y=%s.', ...
        mat2str(size(X)), mat2str(size(Y)));
end

if ~isnan(nx_inferred) && ~isnan(ny_inferred)
    fprintf('Input classified as: %s  (nx=%d, ny=%d)\n', class, nx_inferred, ny_inferred);
else
    fprintf('Input classified as: %s\n', class);
end
end

function tf = y_is_row_constant(Y)
% Structured Cartesian meshgrid has Y constant along each row (rows index Y).
row_std = std(Y, 0, 2);
tf = max(row_std) < 1e-8 * (max(Y(:)) - min(Y(:)));
end
