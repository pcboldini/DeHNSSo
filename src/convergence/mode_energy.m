function [E_tot, E_loc] = mode_energy(StabGrid, Bound, Mode, StabRes, x_loc)
%MODE_ENERGY  Mode-summed perturbation kinetic energy.
%
%   [E_tot, E_loc] = MODE_ENERGY(StabGrid, Bound, Mode, StabRes, x_loc)
%
%   Evaluates the measure DeHNSSo-v1 uses for grid convergence (Eqs. 24 and 25
%   there), summed over the active modes:
%
%       E_tot = int_{x0}^{xb} int_0^H 1/2 (|u'|^2+|v'|^2+|w'|^2) dy dx
%       E_loc = int_0^H       1/2 (|u'|^2+|v'|^2+|w'|^2) dy  at x = x_loc
%
%   The streamwise integral stops at the buffer onset, so the attenuated region
%   never enters the norm. x_loc defaults to 1700, the station Sec. 5.3 and
%   DeHNSSo-v1 both evaluate the local measure at, so the two sections report
%   the same quantity at the same place.
%
%   Unlike the peak amplitude max_y |u'|, this weighs the whole field rather
%   than one wall-normal maximum of one variable, and it weighs the modes by
%   their energy rather than reporting whichever mode differs most.

ny = StabGrid.ny;  nx = StabGrid.nx;
y  = StabGrid.y;  if ~isvector(y); y = y(:,1); end;  y = y(:);
x  = StabGrid.x;  if ~isvector(x); x = x(1,:); end;  x = x(:).';

% The wall-normal grid is stored from the top boundary down to the wall, so it
% has to be flipped before integrating: trapz over a descending coordinate
% returns the integral with the wrong sign.
iy = 1:ny;
if y(end) < y(1); iy = ny:-1:1; y = y(iy); end

ib = min(Bound.ibuf, nx);            % buffer onset (ibuf = nx+1 when off)
xr = 1:ib;

phys = Mode.RunJ(Mode.RunJ >= round(Mode.nf/2));
e_xy = zeros(ny, numel(xr));
for j = phys(:).'
    u = reshape(StabRes.u(j,:,:), ny, nx);
    v = reshape(StabRes.v(j,:,:), ny, nx);
    w = reshape(StabRes.w(j,:,:), ny, nx);
    e_xy = e_xy + 0.5*(abs(u(iy,xr)).^2 + abs(v(iy,xr)).^2 + abs(w(iy,xr)).^2);
end

e_x   = trapz(y, e_xy, 1);           % wall-normal integral, one value per station
E_tot = trapz(x(xr), e_x);

if nargin < 5 || isempty(x_loc); x_loc = 1700; end
assert(x_loc >= x(1) && x_loc <= x(ib), ...
    'x_loc = %g is outside the buffer-free range [%g, %g]', x_loc, x(1), x(ib));
[~, i_loc] = min(abs(x(xr) - x_loc));
E_loc = e_x(i_loc);
end
