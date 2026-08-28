function report_clustering(mug, sig, ag, coord)
%REPORT_CLUSTERING  Print the streamwise clustering summary (all resamplers).
%
%   report_clustering(mug, sig, ag, coord)
%
%   Reports the requested clustering and, when coord is given, what the mesh
%   actually does: the finest dxi against the equidistant one, and the worst
%   ratio between neighbouring dxi (the stretching, 1.00 = none).

if isempty(mug) || all(ag(:) == 1)
    fprintf('xi distribution: equidistant (no refinement)\n');
    return
end

mug = mug(:).';
if isscalar(sig);  sig = sig * ones(size(mug));  end
if isscalar(ag);   ag  = ag  * ones(size(mug));  end

fprintf('xi distribution: Gaussian clustering, %d centre%s\n', ...
        numel(mug), repmat('s', 1, numel(mug) > 1));
for k = 1:numel(mug)
    fprintf('    centre %d: x_c = %-8.3f  sig = %-5.3f  a_g = %-5.3f  (dx ratio = %.2f)\n', ...
            k, mug(k), sig(k), ag(k), 1 / ag(k));
end
% Measured, from the stations actually produced: the finest spacing against the
% equidistant one, and the worst ratio between neighbouring cells.
if nargin >= 4 && ~isempty(coord) && numel(coord) > 2
    d  = diff(coord(:).');
    du = (coord(end) - coord(1)) / (numel(coord) - 1);
    r  = d(2:end) ./ d(1:end-1);
    fprintf('    measured: finest dxi = %.2f x equidistant (%.2fx finer) ; stretching = %.4f\n', ...
            min(d)/du, du/min(d), max(max(r), max(1./r)));
end
end
