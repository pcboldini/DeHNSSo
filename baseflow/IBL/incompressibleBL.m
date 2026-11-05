function [BL] = incompressibleBL(BL)
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
% Non-similar incompressible boundary layer (implicit march).
%
%   BL = incompressibleBL(BL)
%
%   Marches the boundary layer on [S,L] x [0,H] with an implicit second-order
%   scheme in x and Chebyshev collocation in y. Inflow is a Falkner-Skan-Cooke
%   profile at the local Hartree parameter; the nonlinear momentum/continuity
%   system is solved by simple iteration to a 1e-10 threshold on u and v.
%   Continuity is integrated spectrally from the wall, where v = vw.
%   First/second streamwise derivatives predict the next station to speed
%   convergence, and laminar separation is checked at each station (returns if hit).
%
%   Input struct BL: S, L, H, y_i, nu, nx, ny, Ue (1 x nx), We (scalar),
%   optional vw (1 x nx) wall-normal wall velocity, negative for suction.
%   Output adds: u, v, w (ny x nx), y (Chebyshev nodes), derivatives if queried.

%% Inputs
% Name size         unit   explanation
% BL   (1,1) struct [-]    boundary-layer profiles and derivatives

%% Outputs
% Name size         unit   explanation
% BL   (1,1) struct [-]    boundary-layer profiles and derivatives

fprintf('Starting IBL solver.\n');

%% Rename input variables as main structure

nx = BL.nx; ny = BL.ny; S = BL.S; L = BL.L; 
y_i = BL.y_i; H = BL.H; nu = BL.nu; Ue = BL.Ue; We = BL.We;

%% Define convergence error threshold for both u and v

thres = 1e-10;

%% Initialize fields

u = zeros(ny,nx);
v = zeros(ny,nx);
w = zeros(ny,nx);

if Ue(1) == 0
    Ue(1) = 10^-4;          % Incoming (a non-zero value is needed)
end

if ~isscalar(We)
    error(['The input freestream spanwise velocity We is not a scalar, but: [',num2str(size(We)),']']);
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

dx = (L-S)/(nx-1);

[yvec,DM] = chebdif(ny,2);
[y,D1,D2] = MappingMalik(H,y_i,yvec',DM(:,:,1),DM(:,:,2));
clear DM yvec

I = eye(size(D1)); 

%% Initialize inflow with a similar solution at equivalent Hartree parameter

dUe = -FD1d4o(Ue,dx);

m = (dUe(1))*S/Ue(1);       % Falkner-Skan parameter
m = min(m,1);               % Limit m to m<=1 (stagnation)

beta = 2*m/(m+1);           % Hartree parameter
beta = max(beta,-0.19884);  % Limit beta to beta>=-0.19884 (separation)

fprintf(' Using Falkner-Skan inflow profile with beta=%.3f.\n',beta)
eta_inf = 30; % Similar flow freestream condition is invariant to the case
eta_n = 2000;

[y_init,u_init,v_init] = FalknerSkanCookeReal(eta_inf,eta_n,beta,nu,S,Ue(1),true);

u(:,1) = interp1(y_init,u_init,y,'spline',u_init(end))';
v(:,1) = interp1(y_init,v_init,y,'spline','extrap')';

u(:,2) = u(:,1);
v(:,2) = v(:,1);

%% Initialize LHS & RHS including BCs
LHS (1,1)   = 1; % Dirichlet type BC at freestream
LHS (ny,ny) = 1; % Dirichlet type BC at wall
RHSu(ny,1)  = 0; % No-slip u
RHSw(ny,1)  = 0; % No-slip w

% Continuity dv/dy = -du/dx, integrated spectrally from the wall
LHSv          = D1;
LHSv(ny,:)    = 0;
LHSv(ny,ny)   = 1; % Dirichlet type BC at wall (suction/blowing)

%% Start up Inflow Solution (1st order accurate in x)

tic
i = 1;
in = 2:ny-1; 

errv=1;
err=1;
count=0;
while errv>thres || err>thres

    count=count+1;

    %% Momentum
    
    LHS(in,:) = diag(u(in,i+1)/dx)*I(in,:) + diag(v(in,i+1))*D1(in,:) - nu*D2(in,:);
    RHSu(in)  =      u(in,i+1).*u(in,i)/dx + Ue(i+1)*(Ue(i+1)-Ue(i))/dx;
    RHSu(1)  = Ue(i+1); 
    unew     = LHS\RHSu;
    err = max(abs( unew-u(:,i+1) )) / Ue(i+1);
    u(:,i+1) = unew;
    
    %% Continuity
    
    vold=v(:,i+1);
    RHSv     = -(u(:,i+1)-u(:,i))/dx;
    RHSv(ny) = vw(i+1);
    v(:,i+1) = LHSv\RHSv;
    errv = max(abs( vold-v(:,i+1) )) / Ue(i+1);
    
end
fprintf(' Iterations for IC: %u. Remaining error: %3.2g\n',count,err);

% Check for laminar separation
du = (D1*u(:,i+1))/(Ue(i+1)/H);

if du(end)<0
    fprintf('\n Separating inflow detected!\n Ending BL analysis.\n')
    u = u(:,1:2); 
    v = v(:,1:2);

    RHSw(in) = u(in,i+1).*w(in,i)/dx;
    RHSw(1) = We; 
    w(:,i+1) = LHS\RHSw;
    w = w(:,1:2);
    return
end

RHSw(in) = u(in,i+1).*w(in,i)/dx;
RHSw(1) = We; 
w(:,i+1) = LHS\RHSw;

% Estimate next station using finite difference derivatives (first order)
v(:,3) = 2*v(:,2)-v(:,1);
u(:,3) = 2*u(:,2)-u(:,1);

%% March in the x direction  (second order in x)
breaker = false;
percentage = 1;
fprintf(' Starting marching method.')
fprintf(' Progress: %.0f percent',0)
for i = 2:nx-1   
    
    errv=1;
    err=1;
    count=0;

    while errv>thres || err>thres 
        count=count+1;

        %% Momentum
        
        LHS(in,:) = diag(3*u(in,i+1)/(2*dx))*I(in,:) + diag(v(in,i+1))*D1(in,:) - nu*D2(in,:);        
        RHSu(in)  =        u(in,i+1).*(4*u(in,i)-u(in,i-1))/(2*dx) ...
                        + Ue(i+1)*(3*Ue(i+1)-4*Ue(i)+Ue(i-1))/(2*dx);
        RHSu(1) = Ue(i+1);
        unew     = LHS\RHSu;
        err = max(abs( unew-u(:,i+1) )) / Ue(i+1);
        u(:,i+1) = unew; 
        
        %% Continuity
        
        vold=v(:,i+1);
        RHSv     = -(3*u(:,i+1)-4*u(:,i)+u(:,i-1))/(2*dx);
        RHSv(ny) = vw(i+1);
        v(:,i+1) = LHSv\RHSv;
        errv = max(abs( v(:,i+1)-vold )) / Ue(i+1);
    
        %% Check for separation

        du = (D1*u(:,i+1))/(Ue(i+1)/H);
    
        if du(end)<0
            fprintf('\n Separation detected at x/c=%.2f with Ue=%.3f.\n Ending BL analysis.\n',i/nx,Ue(i))
            breaker = true;
            break
        end
    end
    
    RHSw(in) = u(in,i+1).*(4*w(in,i)-w(in,i-1))/(2*dx);
    RHSw(1) = We; 
    w(:,i+1) = LHS\RHSw;
    
    if breaker
        nx=i;
        u = u(:,1:i); 
        v = v(:,1:i);
        w = w(:,1:i);
        break
    end
    
    if 100*(i+1)/nx == percentage
        for j = 1:ceil(log10(percentage))+9; fprintf('\b'); end
        fprintf('%.0f percent ',percentage)
        percentage = percentage + 1;
    end
    
    % Estimate next station usinge finite difference derivatives (second order)
    v(:,i+2) = (5*v(:,i+1)-4*v(:,i)+v(:,i-1))/2;
    u(:,i+2) = (5*u(:,i+1)-4*u(:,i)+u(:,i-1))/2;
    
end

% Remove next station estimation if full solution was obtained
if i==nx-1
    v=v(:,1:end-1);
    u=u(:,1:end-1);
    fprintf('\n')
end

fprintf(' Computation concluded in %.2fs!\n',toc)

%% Outputs

for i=1:nx
    Du(:,i)=D1*u(:,i);
    DDu(:,i)=D2*u(:,i);
    Dv(:,i)=D1*v(:,i);
    DDv(:,i)=D2*v(:,i);
    Dw(:,i)=D1*w(:,i);
    DDw(:,i)=D2*w(:,i);    
end

dudx = -FD1d4o(u, dx);
dvdx = -FD1d4o(v, dx);
dwdx = -FD1d4o(w, dx);

%% Rename output variables as main structure

BL.u = u; BL.v = v; BL.w = w; BL.y = y; BL.x = S + (0:nx-1)*dx;
BL.Du = Du; BL.DDu = DDu; BL.Dv = Dv; BL.DDv = DDv; BL.Dw = Dw; BL.DDw = DDw;
BL.dudx = dudx; BL.dvdx = dvdx; BL.dwdx = dwdx;

end
