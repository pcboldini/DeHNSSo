function f = prepare_forcing(spec)
%PREPARE_FORCING  Wall-forcing spec -> refinement extent + profile.
%
%   f = prepare_forcing(spec)
%
%   Turns a description of what is at the wall into the centre / half-width used
%   to cluster the mesh and a unit-peak profile f.eval(x). spec is either a
%   parameterized actuator (fields xc, w, optional shape = 'cosine') or a
%   distribution (fields x, v), whose xc/w are taken from the support of v.
%
%   Outputs:
%     f.xc, f.w   refinement extent (physical x)
%     f.eval      unit-peak complex profile, 0 outside [xc-w, xc+w]
%     f.kind      'cosine' or 'distribution'

if isfield(spec, 'v') && isfield(spec, 'x')
    xd = spec.x(:).';
    vd = spec.v(:).';
    assert(numel(xd) == numel(vd) && numel(xd) >= 2, ...
        'prepare_forcing: spec.x and spec.v must be equal-length vectors (>= 2).');
    peak = max(abs(vd));
    assert(peak > 0, 'prepare_forcing: spec.v is all zero.');
    act  = find(abs(vd) > 1e-6 * peak);                    % support
    f.xc = 0.5 * (xd(act(1)) + xd(act(end)));
    f.w  = 0.5 * (xd(act(end)) - xd(act(1)));
    assert(f.w > 0, 'prepare_forcing: degenerate support.');
    f.eval = @(x) interp1(xd, vd / peak, x, 'pchip', 0);   % unit peak, phase kept
    f.kind = 'distribution';
    return
end

assert(isfield(spec,'xc') && isfield(spec,'w') && ~isempty(spec.xc) && ~isempty(spec.w), ...
    'prepare_forcing: actuator spec needs xc and w (or x and v for a distribution).');
shape = 'cosine';
if isfield(spec,'shape') && ~isempty(spec.shape); shape = lower(spec.shape); end
assert(strcmp(shape,'cosine'), 'prepare_forcing: unknown shape ''%s'' (only ''cosine'').', shape);

xc = spec.xc;  w = spec.w;
f.xc   = xc;
f.w    = w;
f.eval = @(x) (abs(x - xc) <= w) .* 0.5 .* (1 + cos(pi * (x - xc) / w));
f.kind = shape;
end
