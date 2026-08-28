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
%   lo, hi bound the xi range and n is the station count. mug gives the centres in
%   physical x, [] meaning uniform spacing; sig and ag are the Gaussian width and
%   the peak spacing factor (ag < 1 denser), each broadcast if scalar; x_lo, x_hi
%   are the physical x at lo and hi. Returns 1 x n stations.

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

% Spacing factor fg(t) on the normalised axis t, the strongest reduction from any
% single centre. fg is evaluated on a fine grid in t, NOT on the station index:
% the stations are what we are about to place, so using them as the abscissa puts
% the finest cell away from mug once the spacing varies.
t  = linspace(-1, 1, 20001);
fg = ones(size(t));
for k = 1:numel(mug)
    fg = min(fg, 1 - (1 - ag(k)) * exp(-0.5 * ((t - mug_n(k)) / sig(k)).^2));
end

% Place the stations by inverting the cumulative station density 1/fg, so the
% local spacing follows fg and the finest cell lands on mug.
C     = cumtrapz(t, 1 ./ fg);
C     = (C - C(1)) / (C(end) - C(1));
tn    = interp1(C, t, linspace(0, 1, n), 'linear');
tn(1) = -1;  tn(end) = 1;
coord = lo + (tn + 1) / 2 * (hi - lo);
end
