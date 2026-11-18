function [prof, alpha0, out] = solver_IPSE(omega, beta, StabGrid, Opt)
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
% Incompressible linear PSE (parabolized stability equations) inflow condition.
%
%   [prof, alpha0, out] = solver_IPSE(omega, beta, StabGrid, Opt)
%
%   A self-similar precursor stretches the inflow profile onto [Rstart, R0]
%   (with Opt.pse_precursor = false the stations of StabGrid are used). The
%   local eigenmode of the PSE operator seeds a first-order backward (BDF1)
%   march to R0, with alpha closed at every station by the Chang & Malik
%   (1993) normalization over the velocity components. State [p;u;v;w]. If
%   dx < 1/|alpha_r| the step is stabilized (de Groot & Kotsonis), unless
%   Opt.pse_stab = false. The march has pse_nx stations, by default the
%   pse_npre = 60 precursor stations. The precursor assumes zero pressure
%   gradient: for inflows with a strong pressure gradient (wings) use
%   IC='ILST'. Same structure as solver_CPSE.

%% Inputs
% Name     size unit   explanation
% omega         [-]    angular frequency of the mode
% beta          [-]    spanwise wavenumber of the mode
% StabGrid      [-]    grid/base-flow struct (U, V, W, y, dzdy, d2zdy2, Re, ny)
% Opt           [-]    options (pse_Rstart, pse_npre, pse_nx, pse_maxit, pse_tol,
%                      pse_alpha0, pse_precursor, pse_stab)

%% Outputs
% Name   size unit   explanation
% prof        [-]    inflow shape at R0, normalized by max|u| (u,v,w,p)
% alpha0      [-]    complex streamwise wavenumber at R0
% out         [-]    march history (Rv, x, alpha, phi)

%% Options and defaults
if nargin < 4; Opt = struct(); end
Opt.pse_Rstart = getopt(Opt,'pse_Rstart', max(0.6*StabGrid.Re, StabGrid.Re-200));
Opt.pse_npre   = getopt(Opt,'pse_npre',   60);
Opt.pse_maxit  = getopt(Opt,'pse_maxit',  40);
Opt.pse_tol    = getopt(Opt,'pse_tol',    1e-9);
Opt.pse_stab   = getopt(Opt,'pse_stab',   true);

%% Precursor and march
if getopt(Opt,'pse_precursor',true)
    SG = precursor(StabGrid, Opt.pse_Rstart, Opt.pse_npre);
else
    SG = StabGrid;
end
[prof, alpha0, out] = march(SG, omega, beta, Opt);

end

%% Local functions

function SG = precursor(G, Rstart, Np)
%PRECURSOR  Self-similar upstream extension of the inflow plane.
%   Stretches station 1 of G onto Np stations over [Rstart, R0] (y -> y*R0/R)
%   and recovers V from continuity.

R0 = G.Re;  ny = G.ny;  y = G.y(:);
Rv = linspace(Rstart, R0, Np);
xn = Rv.^2 / R0;

% Edge node, opposite the wall
[~,kw] = min(G.U(:,1));  ed = 1;  if kw == 1; ed = ny; end

SG.U = zeros(ny,Np);  SG.W = zeros(ny,Np);
for i = 1:Np
    ys = y*(R0/Rv(i));
    SG.U(:,i) = interp1(y,G.U(:,1),ys,'pchip',G.U(ed,1));
    SG.W(:,i) = interp1(y,G.W(:,1),ys,'pchip',G.W(ed,1));
end
SG.V = contV(y, xn, SG.U);

SG.ny = ny;  SG.y = y;  SG.Re = R0;  SG.x = xn;
SG.dzdy = G.dzdy;  SG.d2zdy2 = G.d2zdy2;
SG.zpg = true;                                     % zero pressure gradient

end

function V = contV(y, x, U)
%CONTV  Wall-normal velocity from continuity: v = -int_wall^y du/dx dy', v(wall) = 0.

[ny,nx] = size(U);
dUdx = zeros(ny,nx);
for k = 1:ny; dUdx(k,:) = gradient(U(k,:), x); end

[ys,is] = sort(y(:));
V = zeros(ny,nx);
for i = 1:nx
    I = zeros(ny,1);
    I(is) = cumtrapz(ys, dUdx(is,i));
    V(:,i) = -I;
end

end

function [prof, alpha0, out] = march(SG, omega, beta, Opt)
%MARCH  BDF1 PSE march from Rstart to R0 with the alpha closure.

N = SG.ny;  Re = SG.Re;  y = SG.y(:);

% Chebyshev derivatives, mapped to physical y
[~,DM] = chebdif(N,2);  Dn = DM(:,:,1);
D1 = diag(SG.dzdy(:))*Dn;
D2 = diag(SG.d2zdy2(:))*Dn + diag(SG.dzdy(:).^2)*(Dn*Dn);

%% Stations and base flow
R_fine = sqrt(SG.x(:).'*Re);  Rstart = R_fine(1);  R0 = R_fine(end);
nxp = getopt(Opt,'pse_nx', Opt.pse_npre);
if isfield(SG,'zpg') && SG.zpg
    Rv = linspace(Rstart, R0, nxp);                % uniform in R (precursor)
else
    Rv = sqrt(linspace(Rstart^2, R0^2, nxp));      % uniform in x (given stations)
end
xv = Rv.^2 / Re;
bf = interp_bf(SG, R_fine, Rv, D1);

%% Seed: eigenmode of L0 + a*L1 + a^2*L2 (a = i*alpha) at the first station
E0 = operator(bf,1, 0,omega,beta,D1,D2,Re);
Ep = operator(bf,1,+1,omega,beta,D1,D2,Re);
Em = operator(bf,1,-1,omega,beta,D1,D2,Re);
[L0,L1,L2] = eig_bc(E0, (Ep-Em)/2, (Ep+Em)/2-E0, N, D1, D2, Re);
[V,av] = polyeig(L0,L1,L2);
a = -1i*av;  k = isfinite(a) & abs(a) < 1e3;  a = a(k);  V = V(:,k);

% Select the physical mode (alpha_r < 0 kept for stationary crossflow)
Uc = bf.U(:,1);
kd = find(Uc <= 0.99*max(Uc), 1);  if isempty(kd); kd = 1; end
if ~isempty(getopt(Opt,'pse_alpha0',[]))
    [~,i] = min(abs(a - Opt.pse_alpha0));  alpha = a(i);
else
    [alpha,i] = EVfilter(a, V(2*N+1:3*N,:), y, omega, beta, y(kd));
end
phi = V(:,i);

%% Promote the local eigenmode to a PSE eigenmode
dx = xv(2) - xv(1);
for it = 1:20
    [E,G] = operator(bf,1,1i*alpha,omega,beta,D1,D2,Re);
    [A,r] = step_bc(E+G/dx, (G/dx)*phi, N, D1, D2, Re);
    pn = A\r;
    da = closure(pn, (pn-phi)/dx, y, N);
    alpha = alpha + da;  phi = pn;
    if abs(da) < 5e-6; break; end
end

%% March to R0
al  = zeros(1,nxp);           al(1)    = alpha;
PHI = zeros(numel(phi),nxp);  PHI(:,1) = phi;
for j = 2:nxp
    dx = xv(j) - xv(j-1);  pp = PHI(:,j-1);
    if j >= 3; alpha = 2*al(j-1) - al(j-2); else; alpha = al(j-1); end
    for it = 1:Opt.pse_maxit
        [E,G] = operator(bf,j,1i*alpha,omega,beta,D1,D2,Re);
        % Stabilization for steps shorter than the wavelength scale
        s = 1/abs(real(alpha)) - dx;
        if Opt.pse_stab && real(alpha) ~= 0 && s > 0; G = G + s*E; end
        [A,r] = step_bc(E+G/dx, (G/dx)*pp, N, D1, D2, Re);
        pn = A\r;
        da = closure(pn, (pn-pp)/dx, y, N);
        alpha = alpha + da;
        if abs(da) < Opt.pse_tol; break; end
    end
    al(j) = alpha;  PHI(:,j) = pn;
end
out = struct('Rv',Rv, 'x',xv, 'alpha',al, 'phi',PHI);
alpha0 = al(end);

%% Inflow shape at R0, normalized by max|u|
q = PHI(:,end);
p = q(1:N);  u = q(N+1:2*N);  v = q(2*N+1:3*N);  w = q(3*N+1:4*N);
un = max(abs(u));  if un == 0; un = 1; end
prof = struct('u',(u/un).', 'v',(v/un).', 'w',(w/un).', 'p',(p/un).');

end

function bf = interp_bf(SG, R_fine, Rv, D1)
%INTERP_BF  Base flow and its derivatives on the PSE stations.
%   Wall-normal derivatives from SG if given (dyU, ...), else D1; streamwise by FD in x.

xf = R_fine.^2 / SG.Re;  xv = Rv.^2 / SG.Re;
ip  = @(F) interp1(R_fine(:), F.', Rv(:), 'pchip', 'extrap').';
ipx = @(F) interp1(xf(:), ddx(F,xf).', xv(:), 'pchip', 'extrap').';

bf.U  = ip(SG.U);   bf.V  = ip(SG.V);   bf.W  = ip(SG.W);
bf.Ux = ipx(SG.U);  bf.Vx = ipx(SG.V);  bf.Wx = ipx(SG.W);
if isfield(SG,'dyU')
    bf.Uy = ip(SG.dyU);  bf.Vy = ip(SG.dyV);  bf.Wy = ip(SG.dyW);
else
    bf.Uy = D1*bf.U;  bf.Vy = D1*bf.V;  bf.Wy = D1*bf.W;
end

end

function dF = ddx(F, x)
%DDX  Streamwise gradient of each row of F.
dF = zeros(size(F));
for k = 1:size(F,1); dF(k,:) = gradient(F(k,:), x); end
end

function [E,G] = operator(bf, i, a, omega, beta, D1, D2, Re)
%OPERATOR  E and G of E q + G q_x = 0 at station i, a = i*alpha, state [p;u;v;w].
%   Rows: continuity, x-, y-, z-momentum.

N = size(D1,1);  I = eye(N);  Z = zeros(N);  iu = 1i;
U  = diag(bf.U(:,i));   Ux = diag(bf.Ux(:,i));  Uy = diag(bf.Uy(:,i));
V  = diag(bf.V(:,i));   Vx = diag(bf.Vx(:,i));  Vy = diag(bf.Vy(:,i));
W  = diag(bf.W(:,i));   Wx = diag(bf.Wx(:,i));  Wy = diag(bf.Wy(:,i));

% Convection-diffusion
CD = -iu*omega*I + a*U + iu*beta*W + ((beta^2 - a^2)/Re)*I + V*D1 - D2/Re;

E = [Z          a*I      D1       iu*beta*I
     a*I        CD+Ux    Uy       Z
     D1         Vx       CD+Vy    Z
     iu*beta*I  Wx       Wy       CD];

Gm = U - (2*a/Re)*I;
G = [Z  I   Z   Z
     I  Gm  Z   Z
     Z  Z   Gm  Z
     Z  Z   Z   Gm];

end

function [A,r] = step_bc(A, r, N, D1, D2, Re)
%STEP_BC  Boundary rows of the marched step.
%   u,v,w = 0 at the freestream (row 1) and wall (row N); the continuity rows
%   there become -v_yy/Re + p_y = 0.
for blk = 1:3
    for row = [blk*N+1, blk*N+N]
        A(row,:) = 0;  A(row,row) = 1;  r(row) = 0;
    end
end
for row = [1 N]
    A(row,:) = 0;  A(row,1:N) = D1(row,:);  A(row,2*N+1:3*N) = -D2(row,:)/Re;  r(row) = 0;
end
end

function [L0,L1,L2] = eig_bc(L0, L1, L2, N, D1, D2, Re)
%EIG_BC  Boundary rows of the seed eigenproblem, as in step_bc.
for blk = 1:3
    rf = blk*N+1;  rw = blk*N+N;
    L0(rf,:) = 0;  L0(rf,rf) = 1;
    L0(rw,:) = 0;  L0(rw,rw) = 1;
    L1([rf rw],:) = 0;  L2([rf rw],:) = 0;
end
for row = [1 N]
    L0(row,:) = 0;  L0(row,1:N) = D1(row,:);  L0(row,2*N+1:3*N) = -D2(row,:)/Re;
    L1(row,:) = 0;  L2(row,:) = 0;
end
end

function da = closure(phi, dphi, y, N)
%CLOSURE  Chang & Malik normalization, int q^H q_x dy = 0 over u,v,w.
num = 0;  den = 0;
for blk = 1:3
    q  = phi(blk*N+1:blk*N+N);
    qx = dphi(blk*N+1:blk*N+N);
    num = num + trapz(y, conj(q).*qx);
    den = den + trapz(y, abs(q).^2);
end
if abs(den) < 1e-20; da = 0; else; da = -1i*num/den; end
end

function v = getopt(O, f, dv)
%GETOPT  O.f if present and non-empty, else dv.
if isfield(O,f) && ~isempty(O.(f)); v = O.(f); else; v = dv; end
end
