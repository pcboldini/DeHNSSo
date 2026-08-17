function report_clustering(mug, sig, ag)
%REPORT_CLUSTERING  Print the streamwise clustering summary (both resamplers).
%
%   report_clustering(mug, sig, ag)

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
end
