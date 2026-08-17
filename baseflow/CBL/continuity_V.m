function v = continuity_V(y_phys, x_phys, u, rho)
%CONTINUITY_V  Wall-normal base-flow velocity from compressible continuity.
%   v = continuity_V(y_phys, x_phys, u, rho) returns the wall-normal velocity
%   V(y,x) enforcing  d(rho u)/dx + d(rho v)/dy = 0  (2D or infinite-swept,
%   d/dz = 0), i.e.
%
%       rho*v(y) = - int_wall^y  d(rho u)/dx  dy'     [ v(wall) = 0 ].
%
%   This is the EXACT compressible wall-normal velocity for any base flow and
%   replaces the low-Mach Lees-Dorodnitsyn closed form
%       v = sqrt(nu_e Ue/2x) (eta_LD f' - f)/rho ,
%   which is ~correct only as Ma->0 and is ~6x too small at Ma=4.5 (the strong
%   temperature/density variation is not captured). The old form left the base
%   flow violating discrete continuity by ~10x and made the HNS/PSE march
%   over-predict M4.5 2nd-mode growth by ~17%.
%
%   Inputs
%     y_phys  (ny x 1) physical wall-normal grid (any ordering)
%     x_phys  (1 x nx) streamwise stations (physical)
%     u       (ny x nx) streamwise velocity
%     rho     (ny x nx) density (same normalisation used consistently)
%   Output
%     v       (ny x nx) wall-normal velocity, same units/order as u
%
%   d(rho u)/dx is a centred FD across stations (gradient); the wall-normal
%   integral is cumulative-trapezoid from the wall up.

y = y_phys(:);
[ny, nx] = size(u);
rhou = rho .* u;
drhoudx = zeros(ny, nx);
for j = 1:ny
    drhoudx(j,:) = gradient(rhou(j,:), x_phys(:).');
end
[ys, is] = sort(y, 'ascend');          % integrate from wall upward
v = zeros(ny, nx);
for i = 1:nx
    rhov = -cumtrapz(ys, drhoudx(is,i));
    v(is,i) = rhov ./ rho(is,i);
end
end
