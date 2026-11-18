function [prof, alpha0, out] = solver_CPSE(omega, beta, StabGrid, Opt)
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
% Compressible linear PSE (parabolized stability equations) inflow condition (Chang & Malik 1993).
%
%   [prof, alpha0, out] = solver_CPSE(omega, beta, StabGrid, Opt)
%
%   A self-similar precursor stretches the inflow profile onto [Rstart, R0]
%   (with Opt.pse_precursor = false the stations of StabGrid are used, with
%   V and the base pressure gradient from StabGrid). The local eigenmode of
%   the PSE operator seeds a first-order backward (BDF1) march to R0, with
%   alpha closed at every station by the Chang & Malik normalization over the
%   velocity components. Density form, state [rho;u;v;w;T], pressure
%   eliminated through the ideal-gas equation of state. With Opt.pse_stab =
%   true, steps with dx < 1/|alpha_r| are stabilized (de Groot & Kotsonis);
%   off by default, needed if pse_nx is raised.
%   Same structure as solver_IPSE.

%% Inputs
% Name     size unit   explanation
% omega         [-]    angular frequency of the mode
% beta          [-]    spanwise wavenumber of the mode
% StabGrid      [-]    grid/base-flow struct (U, V, W, rho, T, y, dzdy, d2zdy2,
%                      Re, ny, Ma, Pr, gamma, Sstar)
% Opt           [-]    options (pse_Rstart, pse_npre, pse_nx, pse_maxit, pse_tol,
%                      pse_imag_max, pse_alpha0, mode_select, bc_top,
%                      pse_precursor, pse_stab)

%% Outputs
% Name   size unit   explanation
% prof        [-]    inflow shape at R0, normalized by max|u| (u,v,w,rho,T,p)
% alpha0      [-]    complex streamwise wavenumber at R0
% out         [-]    march history (Rv, x, alpha, phi)

%% Options and defaults
if nargin < 4; Opt = struct(); end
Opt.pse_Rstart  = getopt(Opt,'pse_Rstart', max(0.6*StabGrid.Re, StabGrid.Re-200));
Opt.pse_npre    = getopt(Opt,'pse_npre',    60);
Opt.pse_maxit   = getopt(Opt,'pse_maxit',   40);
Opt.pse_tol     = getopt(Opt,'pse_tol',     1e-9);
Opt.mode_select = getopt(Opt,'mode_select', 'auto');
Opt.pse_stab    = getopt(Opt,'pse_stab',    false);

% Top boundary: Dirichlet unless Neumann
bc = getopt(Opt,'bc_top','H_DR');
if iscell(bc)
    bc = 'H_DR';
    if any(strcmpi(getopt(Opt,'bc_top','H_DR'),'H_NM')); bc = 'H_NM'; end
end
Opt.bc_top = bc;

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

SG.U = zeros(ny,Np);  SG.W = zeros(ny,Np);  SG.rho = zeros(ny,Np);  SG.T = zeros(ny,Np);
for i = 1:Np
    ys = y*(R0/Rv(i));
    SG.U(:,i)   = interp1(y,G.U(:,1),ys,'pchip',G.U(ed,1));
    SG.W(:,i)   = interp1(y,G.W(:,1),ys,'pchip',G.W(ed,1));
    SG.rho(:,i) = interp1(y,G.rho(:,1),ys,'pchip',G.rho(ed,1));
    SG.T(:,i)   = interp1(y,G.T(:,1),ys,'pchip',G.T(ed,1));
end
SG.V = contV(y, xn, SG.U, SG.rho);

SG.ny = ny;  SG.y = y;  SG.Re = R0;  SG.x = xn;
SG.dzdy = G.dzdy;  SG.d2zdy2 = G.d2zdy2;
SG.Ma = G.Ma;  SG.Pr = G.Pr;  SG.gamma = G.gamma;  SG.Sstar = G.Sstar;
SG.zpg = true;                                     % zero pressure gradient

end

function V = contV(y, x, U, R)
%CONTV  Wall-normal velocity from continuity: rho*v = -int_wall^y d(rho u)/dx dy', v(wall) = 0.

[ny,nx] = size(U);
dRUdx = zeros(ny,nx);
for k = 1:ny; dRUdx(k,:) = gradient(R(k,:).*U(k,:), x); end

[ys,is] = sort(y(:));
V = zeros(ny,nx);
for i = 1:nx
    I = zeros(ny,1);
    I(is) = cumtrapz(ys, dRUdx(is,i));
    V(:,i) = -I ./ max(R(:,i), eps);
end

end

function [prof, alpha0, out] = march(SG, omega, beta, Opt)
%MARCH  BDF1 PSE march from Rstart to R0 with the alpha closure.

N = SG.ny;  y = SG.y(:);

% Chebyshev derivatives, mapped to physical y
[~,DM] = chebdif(N,2);  Dn = DM(:,:,1);
D1 = diag(SG.dzdy(:))*Dn;
D2 = diag(SG.d2zdy2(:))*Dn + diag(SG.dzdy(:).^2)*(Dn*Dn);

%% Stations and base flow
R_fine = sqrt(SG.x(:).'*SG.Re);  Rstart = R_fine(1);  R0 = R_fine(end);
nxp = getopt(Opt,'pse_nx', max(6, 1+floor((R0^2-Rstart^2)/SG.Re/(1.5/0.04))));
if isfield(SG,'zpg') && SG.zpg
    Rv = linspace(Rstart, R0, nxp);                % uniform in R (precursor)
else
    Rv = sqrt(linspace(Rstart^2, R0^2, nxp));      % uniform in x (given stations)
end
xv = Rv.^2 / SG.Re;
bf = interp_bf(SG, R_fine, Rv, D1, D2);

%% Seed: eigenmode of L0 + a*L1 + a^2*L2 (a = i*alpha) at the first station
E0 = operator(bf,1, 0,omega,beta,D1,D2,SG);
Ep = operator(bf,1,+1,omega,beta,D1,D2,SG);
Em = operator(bf,1,-1,omega,beta,D1,D2,SG);
[L0,L1,L2] = eig_bc(E0, (Ep-Em)/2, (Ep+Em)/2-E0, N, D1, Opt.bc_top);
[V,av] = polyeig(L0,L1,L2);
a = -1i*av;  k = isfinite(a) & abs(a) < 1e3;  a = a(k);  V = V(:,k);

% Physical band (pse_imag_max caps damping; 0.05 keeps near-neutral modes)
imx = getopt(Opt,'pse_imag_max',0.05);
b = real(a) > 0.005 & real(a) < 0.5 & imag(a) > -0.1 & imag(a) < imx;
if any(b); a = a(b);  V = V(:,b); end

% Select the physical mode
Uc = bf.U(:,1);
kd = find(Uc <= 0.99*max(Uc), 1);  if isempty(kd); kd = 1; end
if ~isempty(getopt(Opt,'pse_alpha0',[]))
    [~,i] = min(abs(a - Opt.pse_alpha0));  alpha = a(i);
else
    [alpha,i] = EVfilter(a, V(2*N+1:3*N,:), y, omega, beta, y(kd), SG.Ma, Opt.mode_select);
end
phi = V(:,i);

%% Promote the local eigenmode to a PSE eigenmode
dx = xv(2) - xv(1);
for it = 1:20
    [E,G] = operator(bf,1,1i*alpha,omega,beta,D1,D2,SG);
    [A,r] = step_bc(E+G/dx, (G/dx)*phi, N);
    pn = A\r;
    da = closure(pn, (pn-phi)/dx, bf.rho(:,1), y, N);
    alpha = alpha + da;  phi = pn;
    if abs(da) < 5e-6; break; end
end

%% March to R0
al  = zeros(1,nxp);   al(1)    = alpha;
PHI = zeros(5*N,nxp); PHI(:,1) = phi;
for j = 2:nxp
    dx = xv(j) - xv(j-1);  pp = PHI(:,j-1);
    if j >= 3; alpha = 2*al(j-1) - al(j-2); else; alpha = al(j-1); end
    for it = 1:Opt.pse_maxit
        [E,G] = operator(bf,j,1i*alpha,omega,beta,D1,D2,SG);
        % Stabilization for steps shorter than the wavelength scale
        s = 1/abs(real(alpha)) - dx;
        if Opt.pse_stab && real(alpha) ~= 0 && s > 0; G = G + s*E; end
        [A,r] = step_bc(E+G/dx, (G/dx)*pp, N);
        pn = A\r;
        da = closure(pn, (pn-pp)/dx, bf.rho(:,j), y, N);
        alpha = alpha + da;
        if abs(da) < Opt.pse_tol; break; end
    end
    al(j) = alpha;  PHI(:,j) = pn;
end
out = struct('Rv',Rv, 'x',xv, 'alpha',al, 'phi',PHI);
alpha0 = al(end);

%% Inflow shape at R0, normalized by max|u|
q  = PHI(:,end);
rh = q(1:N);  u = q(N+1:2*N);  v = q(2*N+1:3*N);  w = q(3*N+1:4*N);  T = q(4*N+1:5*N);
p  = (rh.*bf.T(:,end) + bf.rho(:,end).*T) / (SG.gamma*SG.Ma^2);
un = max(abs(u));  if un == 0; un = 1; end
prof = struct('u',(u/un).', 'v',(v/un).', 'w',(w/un).', ...
              'rho',(rh/un).', 'T',(T/un).', 'p',(p/un).');

end

function bf = interp_bf(SG, R_fine, Rv, D1, D2)
%INTERP_BF  Base flow and its derivatives on the PSE stations.
%   Wall-normal derivatives from SG if given (dyU, ...), else D1/D2; streamwise by FD in x.

xf = R_fine.^2 / SG.Re;  xv = Rv.^2 / SG.Re;
ip  = @(F) interp1(R_fine(:), F.', Rv(:), 'pchip', 'extrap').';
ipx = @(F) interp1(xf(:), ddx(F,xf).', xv(:), 'pchip', 'extrap').';

bf.U  = ip(SG.U);   bf.V  = ip(SG.V);   bf.T   = ip(SG.T);   bf.rho = ip(SG.rho);
bf.Ux = ipx(SG.U);  bf.Vx = ipx(SG.V);  bf.Tx  = ipx(SG.T);  bf.Rx  = ipx(SG.rho);
if isfield(SG,'dyU')
    bf.Uy = ip(SG.dyU);  bf.Uyy = ip(SG.ddyU);  bf.Vy = ip(SG.dyV);
    bf.Ty = ip(SG.dyT);  bf.Tyy = ip(SG.ddyT);  bf.Ry = ip(SG.dyrho);
else
    bf.Uy = D1*bf.U;  bf.Uyy = D2*bf.U;  bf.Vy = D1*bf.V;
    bf.Ty = D1*bf.T;  bf.Tyy = D2*bf.T;  bf.Ry = D1*bf.rho;
end
if isfield(SG,'zpg') && SG.zpg
    bf.Px = zeros(size(bf.U));                     % zero pressure gradient
else
    bf.Px = ipx(SG.rho.*SG.T)/(SG.gamma*SG.Ma^2);  % p = rho T/(gamma Ma^2)
end
[bf.mu, bf.muT, bf.muTT] = sutherland(bf.T, SG.Sstar);

end

function dF = ddx(F, x)
%DDX  Streamwise gradient of each row of F.
dF = zeros(size(F));
for k = 1:size(F,1); dF(k,:) = gradient(F(k,:), x); end
end

function [E,G] = operator(bf, i, a, omega, beta, D1, D2, SG)
%OPERATOR  E and G of E q + G q_x = 0 at station i, a = i*alpha, state [rho;u;v;w;T].
%   Density form of Chang & Malik (1993): calorically perfect gas, Sutherland
%   viscosity, constant Prandtl number (k = mu).

N = size(D1,1);  I = eye(N);  Z = zeros(N);  iu = 1i;
Re = SG.Re;  kM2 = SG.gamma*SG.Ma^2;  Ec = (SG.gamma-1)*SG.Ma^2;  sig = SG.Pr;
dU = diag(bf.U(:,i));  dUy = diag(bf.Uy(:,i));  dUyy = diag(bf.Uyy(:,i));  dUx = diag(bf.Ux(:,i));
dV = diag(bf.V(:,i));  dVy = diag(bf.Vy(:,i));  dVx = diag(bf.Vx(:,i));
dT = diag(bf.T(:,i));  dTy = diag(bf.Ty(:,i));  dTyy = diag(bf.Tyy(:,i));  dTx = diag(bf.Tx(:,i));
dR = diag(bf.rho(:,i));  dRy = diag(bf.Ry(:,i));  dRx = diag(bf.Rx(:,i));
dmu = diag(bf.mu(:,i));  dmuT = diag(bf.muT(:,i));  dmuTT = diag(bf.muTT(:,i));
dmuY = diag(D1*bf.mu(:,i));  dPx = diag(bf.Px(:,i));
dkap = dmu;  dkapT = dmuT;  dkapTT = dmuTT;  dcp = I;  dcpT = Z;  dOmega = I;

E11 = -iu*omega*I+dU*a+dUx+dV*D1+dVy;           E12 = dRx+dR*a; E13 = dRy+dR*D1; E14 = iu*beta*dR; E15 = Z;
E21 = dU*dUx+dV*dUy+(1/kM2)*(dT*a+dTx);
E22 = dR*(-iu*omega*I+dU*a+dV*D1+dUx)-(dmu/Re)*((4/3)*a^2*I-beta^2*I+D2)-(1/Re)*dmuY*D1;
E23 = dR*dUy-(a/Re)*(dmuY+(dmu/3)*D1);           E24 = -iu*beta*dmu*a/(3*Re);
E25 = -(dmuT/Re)*(dUyy+dUy*D1)+(1/kM2)*(dRx+dR*a)-(dmuTT/Re)*dTy*dUy;
E31 = dU*dVx+dV*dVy+(1/kM2)*(dTy+dT*D1);         E32 = -(a/(3*Re))*(dmu*D1-2*dmuY);
E33 = dR*(-iu*omega*I+dU*a+dV*D1+dVy)-(dmu/Re)*(a^2*I-beta^2*I+(4/3)*D2)-(4/(3*Re))*dmuY*D1;
E34 = -(iu*beta/(3*Re))*(dmu*D1-2*dmuY);         E35 = (1/kM2)*(dR*D1+dRy)-(dmuT*a/Re)*dUy;
E41 = iu*beta*dT/kM2;                            E42 = -iu*beta*a*dmu/(3*Re);
E43 = -(iu*beta/Re)*((dmu/3)*D1+dmuY);
E44 = dR*(-iu*omega*I+dU*a+dV*D1)-(dmu/Re)*(a^2*I-(4/3)*beta^2*I+D2)-(1/Re)*dmuY*D1;
E45 = iu*beta*dR/kM2;
E51 = (dcp-(Ec/kM2)*I)*(dU*dTx+dV*dTy)-(Ec*dT/kM2)*(-iu*omega*I+dU*a+dV*D1);
E52 = dcp*dR*dTx-2*Ec*(dmu/Re)*dUy*D1-Ec*dPx;    E53 = dcp*dR*dTy-2*Ec*(dmuT/Re)*dUy*a; E54 = Z;
E55 = dcp*dR*(-iu*omega*I+dU*a+dV*D1)-Ec*(dmuT/Re)*dUy.^2+dR*dcpT*(dU*dTx+dV*dTy) ...
     -(Ec/kM2)*(-iu*omega*dR+dR*dU*a+dR*dV*D1+dV*dRy+dU*dRx) ...
     -(1/(Re*sig))*(dkap*(a^2*I-beta^2*I+D2)+dkapT*(dTyy+2*dTy*D1)+dkapTT*dTy.^2);
E = [E11 E12 E13 E14 E15; E21 E22 E23 E24 E25; E31 E32 E33 E34 E35; E41 E42 E43 E44 E45; E51 E52 E53 E54 E55];

G11 = dU; G12 = dR; G13 = Z; G14 = Z; G15 = Z;
G21 = dOmega*(dT/kM2); G22 = dR*dU-(8/3)*(dmu/Re)*a; G23 = -(1/Re)*(dmuY+(dmu/3)*D1);
G24 = -(iu*beta/(3*Re))*dmu; G25 = dOmega*(dR/kM2);
G31 = Z; G32 = -(1/(3*Re))*(dmu*D1-2*dmuY); G33 = dR*dU-2*(dmu/Re)*a; G34 = Z; G35 = Z;
G41 = Z; G42 = -(iu*beta/(3*Re))*dmu; G43 = Z; G44 = dR*dU-2*(dmu/Re)*a; G45 = Z;
G51 = -Ec*dT*dU/kM2; G52 = Z; G53 = Z; G54 = Z; G55 = dR*dU*(dcp-(Ec/kM2)*I)-2*(dkap/(Re*sig))*a;
G = [G11 G12 G13 G14 G15; G21 G22 G23 G24 G25; G31 G32 G33 G34 G35; G41 G42 G43 G44 G45; G51 G52 G53 G54 G55];

end

function [A,r] = step_bc(A, r, N)
%STEP_BC  Boundary rows of the marched step: u,v,w,T = 0 at the freestream
%   (row 1) and wall (row N); rho free.
for blk = 1:4
    for row = [blk*N+1, blk*N+N]
        A(row,:) = 0;  A(row,row) = 1;  r(row) = 0;
    end
end
end

function [L0,L1,L2] = eig_bc(L0, L1, L2, N, D1, bc)
%EIG_BC  Boundary rows of the seed eigenproblem: Dirichlet u,v,w,T, or Neumann
%   at the freestream for H_NM; rho free.
for blk = 1:4
    rf = blk*N+1;  rw = blk*N+N;
    L0(rf,:) = 0;
    if strcmpi(bc,'H_NM'); L0(rf,blk*N+1:blk*N+N) = D1(1,:); else; L0(rf,rf) = 1; end
    L0(rw,:) = 0;  L0(rw,rw) = 1;
    L1([rf rw],:) = 0;  L2([rf rw],:) = 0;
end
end

function da = closure(phi, dphi, rho, y, N)
%CLOSURE  Chang & Malik normalization, int rho q^H q_x dy = 0 over u,v,w.
w = rho(:);  num = 0;  den = 0;
for blk = 1:3
    q  = phi(blk*N+1:blk*N+N);
    qx = dphi(blk*N+1:blk*N+N);
    num = num + trapz(y, w.*conj(q).*qx);
    den = den + trapz(y, w.*abs(q).^2);
end
if abs(den) < 1e-20; da = 0; else; da = -1i*num/den; end
end

function v = getopt(O, f, dv)
%GETOPT  O.f if present and non-empty, else dv.
if isfield(O,f) && ~isempty(O.(f)); v = O.(f); else; v = dv; end
end
