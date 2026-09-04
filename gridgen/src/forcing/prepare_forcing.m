function f = prepare_forcing(spec)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Wall-forcing spec -> refinement extent + profile.
%
%   f = prepare_forcing(spec)
%
%   Turns a description of what is at the wall into the centre / half-width used
%   to cluster the mesh and a unit-peak profile f.eval(x). spec is either a
%   parameterized actuator (fields xc, w, optional shape = 'cosine') or a
%   distribution (fields x, v), whose xc/w are taken from the support of v.
%
%   Returns f.xc and f.w, the centre and half-width in physical x; f.eval, a
%   unit-peak complex profile that is zero outside [xc-w, xc+w]; and f.kind,
%   'cosine' or 'distribution'.

%% Inputs
% Name size         unit   explanation
% spec (1,1) struct [-]    description of what sits at the wall for one forcing source

%% Outputs
% Name size         unit   explanation
% f    (1,1) struct [-]    refinement extent and wall profile: centre, half-width and shape

% Distribution: the support is where |v| rises above 1e-6 of its peak, and the
% centre and half-width follow from its two ends.
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

% Actuator: raised cosine, unit peak at xc and zero at xc +- w, so w is the
% HALF-width and the footprint spans 2w.
xc = spec.xc;  w = spec.w;
f.xc   = xc;
f.w    = w;
f.eval = @(x) (abs(x - xc) <= w) .* 0.5 .* (1 + cos(pi * (x - xc) / w));
f.kind = shape;
end
