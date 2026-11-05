function [BL] = compressibleBL(BL)
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
% Non-similar compressible boundary layer (implicit march).
%
%   BL = compressibleBL(BL)
%
%   Compressible version of incompressibleBL. Marches the boundary layer on
%   [S,L] x [0,H] with an implicit second-order scheme in x and Chebyshev
%   collocation in y. Perfect gas, Sutherland viscosity, constant pressure
%   across the layer. Inflow is a compressible Falkner-Skan-Cooke profile at the
%   local Hartree parameter; momentum, energy and continuity are solved by simple
%   iteration to a 1e-10 threshold on u, v and T. Continuity is integrated
%   spectrally from the wall, where rho*v = rho_w*vw. The edge temperature and
%   pressure follow from the freestream stagnation state and Ue (isentropic).
%
%   Input struct BL: S, L, H, y_i, nx, ny, Ue (1 x nx), We (scalar),
%   Tinf [K], pinf [Pa], Uinf [m/s] (freestream),
%   optional gamma (1.4), R (287.0587), Pr (0.72), mu_ref (1.716e-5),
%   T_ref (273.15), S_suth (110.4), wall_bc ('adiabatic' or 'isothermal'),
%   Tw [K] (isothermal wall), vw (1 x nx) wall-normal wall velocity, negative for suction.
%   Output adds: u, v, w, T, rho, mu (ny x nx), y (Chebyshev nodes), Te, pe, rhoe,
%   mue (1 x nx), derivatives.

%% Inputs
% Name size         unit   explanation
% BL   (1,1) struct [-]    boundary-layer profiles and derivatives

%% Outputs
% Name size         unit   explanation
% BL   (1,1) struct [-]    boundary-layer profiles and derivatives

fprintf('Starting CBL solver.\n');

%% Rename input variables as main structure

nx = BL.nx; ny = BL.ny; S = BL.S; L = BL.L;
y_i = BL.y_i; H = BL.H; Ue = BL.Ue; We = BL.We;

% Gas properties (defaults: air)
if ~isfield(BL,'gamma');   BL.gamma   = 1.4;         end
if ~isfield(BL,'R');       BL.R       = 287.0587;    end
if ~isfield(BL,'Pr');      BL.Pr      = 0.72;        end
if ~isfield(BL,'mu_ref');  BL.mu_ref  = 1.716e-5;    end
if ~isfield(BL,'T_ref');   BL.T_ref   = 273.15;      end
if ~isfield(BL,'S_suth');  BL.S_suth  = 110.4;       end
if ~isfield(BL,'wall_bc'); BL.wall_bc = 'adiabatic'; end

gamma = BL.gamma; Rgas = BL.R; Pr = BL.Pr;
cp = gamma*Rgas/(gamma-1);

% Sutherland's law
mu_ref = BL.mu_ref; T_ref = BL.T_ref; S_suth = BL.S_suth;
Sutherland = @(T) mu_ref*(T/T_ref).^(3/2)*(T_ref+S_suth)./(T+S_suth);

%% Define convergence error threshold for u, v and T

thres = 1e-10;

%% Initialize fields

u   = zeros(ny,nx);
v   = zeros(ny,nx);
w   = zeros(ny,nx);
T   = zeros(ny,nx);
rho = zeros(ny,nx);

if ~isscalar(We)
    error(['The input freestream spanwise velocity We is not a scalar, but: [',num2str(size(We)),']']);
end

isothermal = strcmpi(BL.wall_bc,'isothermal');
if isothermal
    Tw = BL.Tw.*ones(1,nx);
end

% Wall suction/blowing (zero by default)
if isfield(BL,'vw') && ~isempty(BL.vw)
    vw = BL.vw;
    if numel(vw) ~= nx
        error('BL.vw must have nx = %d entries, not %d.', nx, numel(vw));
    end
    if vw(1) ~= 0 || vw(2) ~= 0
        warning('Suction at the first two stations: the similarity inflow has an impermeable wall.');
    end
else
    vw = zeros(1,nx);
end

%% Edge state from the freestream stagnation state (isentropic)

T0   = BL.Tinf + BL.Uinf^2/(2*cp);
p0   = BL.pinf*(T0/BL.Tinf)^(gamma/(gamma-1));
Te   = T0 - (Ue.^2+We^2)/(2*cp);
pe   = p0*(Te/T0).^(gamma/(gamma-1));
rhoe = pe./(Rgas*Te);
mue  = Sutherland(Te);

dx = (L-S)/(nx-1);

[yvec,DM] = chebdif(ny,2);
[y,D1,D2] = MappingMalik(H,y_i,yvec',DM(:,:,1),DM(:,:,2));
y = y(:);
clear DM yvec

I = eye(size(D1));

%% Initialize inflow with a similar solution at equivalent Hartree parameter

dUe = -FD1d4o(Ue,dx);

m = (dUe(1))*S/Ue(1);       % Falkner-Skan parameter
m = min(m,1);               % Limit m to m<=1 (stagnation)

beta = 2*m/(m+1);           % Hartree parameter
beta = max(beta,-0.19884);  % Limit beta to beta>=-0.19884 (separation)
m = beta/(2-beta);

fprintf(' Using compressible Falkner-Skan inflow profile with beta=%.3f.\n',beta)

FSC.Ma       = Ue(1)/sqrt(gamma*Rgas*Te(1));
FSC.Pr       = Pr;
FSC.gamma    = gamma;
FSC.Sstar    = S_suth/Te(1);
FSC.beta_H   = beta;
FSC.We_ratio = We/Ue(1);
FSC.wall_bc  = BL.wall_bc;
FSC.eta_inf  = 30;
FSC.npts     = 801;
if isothermal
    FSC.Tw = Tw(1)/Te(1);
end
[~,FSC] = compressibleFalknerSkanCooke(FSC);

% compressibleFalknerSkanCooke similarity variable: eta = y sqrt((m+1) Ue/(nu_e x))
eta = y*sqrt((m+1)*Ue(1)*rhoe(1)/(mue(1)*S));
u(:,1)   = Ue(1)*interp1(FSC.eta,FSC.fp,eta,'spline',FSC.fp(end));
w(:,1)   = We*interp1(FSC.eta,FSC.g,eta,'spline',FSC.g(end));
T(:,1)   = Te(1)*interp1(FSC.eta,FSC.theta,eta,'spline',FSC.theta(end));
rho(:,1) = pe(1)./(Rgas*T(:,1));

u(:,2) = u(:,1);
w(:,2) = w(:,1);
T(:,2) = T(:,1);

%% Initialize LHS & RHS including BCs
LHS (1,1)   = 1; % Dirichlet type BC at freestream
LHS (ny,ny) = 1; % Dirichlet type BC at wall
RHSu(ny,1)  = 0; % No-slip u
RHSw(ny,1)  = 0; % No-slip w

LHST      = zeros(ny);
LHST(1,1) = 1;   % Dirichlet type BC at freestream (T = Te)
if isothermal
    LHST(ny,ny) = 1;          % Dirichlet type BC at wall (T = Tw)
else
    LHST(ny,:)  = D1(ny,:);   % Neumann type BC at wall (adiabatic)
end
RHST(ny,1) = 0;

% Continuity d(rho v)/dy = -d(rho u)/dx, integrated spectrally from the wall
LHSv          = D1;
LHSv(ny,:)    = 0;
LHSv(ny,ny)   = 1; % Dirichlet type BC at wall (suction/blowing)

%% Start up Inflow Solution (1st order accurate in x)

tic
i = 1;
in = 2:ny-1;

% Pressure gradient -dpe/dx
dpdx = rhoe(i+1)*Ue(i+1)*(Ue(i+1)-Ue(i))/dx;

errv=1;
err=1;
errT=1;
count=0;
while err>thres || errv>thres || errT>thres

    count=count+1;

    % Gas properties at the current temperature
    rho(:,i+1) = pe(i+1)./(Rgas*T(:,i+1));
    mu = Sutherland(T(:,i+1));
    k  = cp*mu/Pr;

    %% Momentum

    LHS(in,:) = diag(rho(in,i+1).*u(in,i+1)/dx)*I(in,:) + diag(rho(in,i+1).*v(in,i+1))*D1(in,:) ...
              - diag(mu(in))*D2(in,:) - diag(D1(in,:)*mu)*D1(in,:);
    RHSu(in)  = rho(in,i+1).*u(in,i+1).*u(in,i)/dx + dpdx;
    RHSu(1)  = Ue(i+1);
    unew     = LHS\RHSu;
    err = max(abs( unew-u(:,i+1) )) / Ue(i+1);
    u(:,i+1) = unew;

    RHSw(in) = rho(in,i+1).*u(in,i+1).*w(in,i)/dx;
    RHSw(1)  = We;
    w(:,i+1) = LHS\RHSw;

    %% Energy

    LHST(in,:) = cp*diag(rho(in,i+1).*u(in,i+1)/dx)*I(in,:) + cp*diag(rho(in,i+1).*v(in,i+1))*D1(in,:) ...
               - diag(k(in))*D2(in,:) - diag(D1(in,:)*k)*D1(in,:);
    % Convection, pressure work and viscous dissipation
    RHST(in)  = cp*rho(in,i+1).*u(in,i+1).*T(in,i)/dx - u(in,i+1)*dpdx ...
              + mu(in).*((D1(in,:)*u(:,i+1)).^2 + (D1(in,:)*w(:,i+1)).^2);
    RHST(1) = Te(i+1);
    if isothermal
        RHST(ny) = Tw(i+1);
    end
    Tnew = LHST\RHST;
    errT = max(abs( Tnew-T(:,i+1) )) / Te(i+1);
    T(:,i+1)   = Tnew;
    rho(:,i+1) = pe(i+1)./(Rgas*T(:,i+1));

    %% Continuity

    vold=v(:,i+1);
    RHSv     = -(rho(:,i+1).*u(:,i+1)-rho(:,i).*u(:,i))/dx;
    RHSv(ny) = rho(ny,i+1)*vw(i+1);
    v(:,i+1) = (LHSv\RHSv)./rho(:,i+1);
    errv = max(abs( vold-v(:,i+1) )) / Ue(i+1);

end
fprintf(' Iterations for IC: %u. Remaining error: %3.2g\n',count,err);

% Check for laminar separation
du = (D1*u(:,i+1))/(Ue(i+1)/H);

if du(end)<0
    fprintf('\n Separating inflow detected!\n Ending BL analysis.\n')
    u = u(:,1:2); v = v(:,1:2); w = w(:,1:2); T = T(:,1:2); rho = rho(:,1:2);
    return
end

% Estimate next station using finite difference derivatives (first order)
u(:,3) = 2*u(:,2)-u(:,1);
w(:,3) = 2*w(:,2)-w(:,1);
T(:,3) = 2*T(:,2)-T(:,1);
v(:,3) = v(:,2);

%% March in the x direction  (second order in x)
breaker = false;
percentage = 1;
fprintf(' Starting marching method.')
fprintf(' Progress: %.0f percent',0)
for i = 2:nx-1

    % Pressure gradient -dpe/dx
    dpdx = rhoe(i+1)*Ue(i+1)*(3*Ue(i+1)-4*Ue(i)+Ue(i-1))/(2*dx);

    errv=1;
    err=1;
    errT=1;
    count=0;

    while err>thres || errv>thres || errT>thres
        count=count+1;

        % Gas properties at the current temperature
        rho(:,i+1) = pe(i+1)./(Rgas*T(:,i+1));
        mu = Sutherland(T(:,i+1));
        k  = cp*mu/Pr;

        %% Momentum

        LHS(in,:) = diag(3*rho(in,i+1).*u(in,i+1)/(2*dx))*I(in,:) + diag(rho(in,i+1).*v(in,i+1))*D1(in,:) ...
                  - diag(mu(in))*D2(in,:) - diag(D1(in,:)*mu)*D1(in,:);
        RHSu(in)  = rho(in,i+1).*u(in,i+1).*(4*u(in,i)-u(in,i-1))/(2*dx) + dpdx;
        RHSu(1) = Ue(i+1);
        unew     = LHS\RHSu;
        err = max(abs( unew-u(:,i+1) )) / Ue(i+1);
        u(:,i+1) = unew;

        RHSw(in) = rho(in,i+1).*u(in,i+1).*(4*w(in,i)-w(in,i-1))/(2*dx);
        RHSw(1)  = We;
        w(:,i+1) = LHS\RHSw;

        %% Energy

        LHST(in,:) = cp*diag(3*rho(in,i+1).*u(in,i+1)/(2*dx))*I(in,:) + cp*diag(rho(in,i+1).*v(in,i+1))*D1(in,:) ...
                   - diag(k(in))*D2(in,:) - diag(D1(in,:)*k)*D1(in,:);
        % Convection, pressure work and viscous dissipation
        RHST(in)  = cp*rho(in,i+1).*u(in,i+1).*(4*T(in,i)-T(in,i-1))/(2*dx) - u(in,i+1)*dpdx ...
                  + mu(in).*((D1(in,:)*u(:,i+1)).^2 + (D1(in,:)*w(:,i+1)).^2);
        RHST(1) = Te(i+1);
        if isothermal
            RHST(ny) = Tw(i+1);
        end
        Tnew = LHST\RHST;
        errT = max(abs( Tnew-T(:,i+1) )) / Te(i+1);
        T(:,i+1)   = Tnew;
        rho(:,i+1) = pe(i+1)./(Rgas*T(:,i+1));

        %% Continuity

        vold=v(:,i+1);
        RHSv     = -(3*rho(:,i+1).*u(:,i+1)-4*rho(:,i).*u(:,i)+rho(:,i-1).*u(:,i-1))/(2*dx);
        RHSv(ny) = rho(ny,i+1)*vw(i+1);
        v(:,i+1) = (LHSv\RHSv)./rho(:,i+1);
        errv = max(abs( v(:,i+1)-vold )) / Ue(i+1);

        %% Check for separation

        du = (D1*u(:,i+1))/(Ue(i+1)/H);

        if du(end)<0
            fprintf('\n Separation detected at x=%.4f m with Ue=%.3f.\n Ending BL analysis.\n',S+i*dx,Ue(i))
            breaker = true;
            break
        end
    end

    if breaker
        nx=i;
        u = u(:,1:i); v = v(:,1:i); w = w(:,1:i); T = T(:,1:i); rho = rho(:,1:i);
        Te = Te(1:i); pe = pe(1:i); rhoe = rhoe(1:i); mue = mue(1:i);
        break
    end

    if 100*(i+1)/nx == percentage
        for j = 1:ceil(log10(percentage))+9; fprintf('\b'); end
        fprintf('%.0f percent ',percentage)
        percentage = percentage + 1;
    end

    % Estimate next station using finite difference derivatives (second order)
    u(:,i+2) = (5*u(:,i+1)-4*u(:,i)+u(:,i-1))/2;
    v(:,i+2) = (5*v(:,i+1)-4*v(:,i)+v(:,i-1))/2;
    w(:,i+2) = (5*w(:,i+1)-4*w(:,i)+w(:,i-1))/2;
    T(:,i+2) = (5*T(:,i+1)-4*T(:,i)+T(:,i-1))/2;

end

% Remove next station estimation if full solution was obtained
if i==nx-1
    u=u(:,1:end-1);
    v=v(:,1:end-1);
    w=w(:,1:end-1);
    T=T(:,1:end-1);
    rho=rho(:,1:nx);
    fprintf('\n')
end

% Inflow v by extrapolation (not given by the similarity solution)
v(:,1) = 2*v(:,2)-v(:,3);

fprintf(' Computation concluded in %.2fs!\n',toc)

%% Outputs

for i=1:nx
    Du(:,i)=D1*u(:,i);
    DDu(:,i)=D2*u(:,i);
    Dv(:,i)=D1*v(:,i);
    DDv(:,i)=D2*v(:,i);
    Dw(:,i)=D1*w(:,i);
    DDw(:,i)=D2*w(:,i);
    DT(:,i)=D1*T(:,i);
    DDT(:,i)=D2*T(:,i);
end

dudx = -FD1d4o(u, dx);
dvdx = -FD1d4o(v, dx);
dwdx = -FD1d4o(w, dx);
dTdx = -FD1d4o(T, dx);

%% Rename output variables as main structure

BL.u = u; BL.v = v; BL.w = w; BL.T = T; BL.rho = rho; BL.mu = Sutherland(T);
BL.y = y; BL.x = S + (0:nx-1)*dx;
BL.Te = Te; BL.pe = pe; BL.rhoe = rhoe; BL.mue = mue;
BL.Du = Du; BL.DDu = DDu; BL.Dv = Dv; BL.DDv = DDv; BL.Dw = Dw; BL.DDw = DDw;
BL.DT = DT; BL.DDT = DDT;
BL.dudx = dudx; BL.dvdx = dvdx; BL.dwdx = dwdx; BL.dTdx = dTdx;

end
