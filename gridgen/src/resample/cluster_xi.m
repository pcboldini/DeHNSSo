function coord = cluster_xi(lo, hi, n, mug, sig, ag, x_lo, x_hi)
%CLUSTER_XI  Gaussian-clustered stations on [lo, hi].
%
%   coord = cluster_xi(lo, hi, n, mug, sig, ag, x_lo, x_hi)
%
%   Returns n monotonic stations on [lo, hi], the local spacing reduced around
%   one or more centres. Shared by the flat-plate (xi = x) and curvilinear
%   (xi = arc length) resamplers; centres mug are physical x, mapped onto the
%   normalised axis by x_lo -> -1, x_hi -> +1.
%
%   Inputs:
%     lo, hi      — xi range
%     n           — number of stations
%     mug         — centre(s), physical x ([] gives uniform spacing)
%     sig, ag     — Gaussian width / peak spacing factor (ag < 1 denser); scalar broadcasts
%     x_lo, x_hi  — physical x at lo / hi
%   Outputs:
%     coord       — 1 x n stations

if isempty(mug) || all(ag(:) == 1)
    coord = linspace(lo, hi, n);
    return
end

mug = mug(:).';
if isscalar(sig);  sig = sig * ones(size(mug));  end
if isscalar(ag);   ag  = ag  * ones(size(mug));  end
assert(numel(sig) == numel(mug) && numel(ag) == numel(mug), ...
    'cluster_xi: mug/sig/ag must be scalars or same-length vectors.');

mug_n = 2 * (mug - x_lo) / (x_hi - x_lo) - 1;          % centres on [-1, 1]

% Spacing factor fg(xi): the strongest reduction from any single centre.
XC = linspace(-1, 1, n);
fg = ones(size(XC));
for k = 1:numel(mug)
    fg = min(fg, 1 - (1 - ag(k)) * exp(-0.5 * ((XC - mug_n(k)) / sig(k)).^2));
end

% Accumulate into stations and normalise to [lo, hi].
c     = [0, cumsum(fg(2:end))];
coord = lo + (c / c(end)) * (hi - lo);
end
