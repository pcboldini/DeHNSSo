function [y,u,v] = FalknerSkanCookeReal(eta_inf,n,beta,nu,x,Ue,halfsine,tolerance)
%FALKNERSKANCOOKEREAL  Falkner-Skan-Cooke solution in physical coordinates.
%
%   [y,u,v] = FalknerSkanCookeReal(eta_inf, n, beta, nu, x, Ue, halfsine, tolerance)
%
%   Solves the similarity problem (FalknerSkanCooke) and maps it to physical
%   coordinates at station x, with m = beta/(2-beta):
%     y = eta sqrt(2 nu x / ((m+1) Ue)),   u = f' Ue,
%     v = -sqrt(nu Ue/(2 x (m+1))) ((m+1) f + (m-1) f' eta).
%
%   Inputs:
%     eta_inf, n, halfsine, tolerance  — passed to FalknerSkanCooke
%     beta  — Hartree parameter,  nu — kinematic viscosity
%     x     — station,            Ue — local streamwise edge velocity
%   Outputs:
%     y, u, v  — physical wall-normal grid and velocity components


if nargin < 8
    halfsine = false;
end
if nargin < 9
    tolerance = 1e-8;
end
if beta==2
    error('Beta cant be equal to 2!')
end

% Get similar solution
[sol] = FalknerSkanCooke(eta_inf,n,beta,halfsine,tolerance);

% Transform Hartree to Falkner-Skan parameter
m = beta/(2-beta); % beta = 2*m/(1+m)

% Scale to real coordinates
y = sol.x * sqrt(nu*x/Ue * 2/(m+1));
u = sol.y(2,:) * Ue;
v = -sqrt(nu*Ue/(2*x*(m+1))) * ((m+1)*sol.y(1,:)+(m-1)*sol.y(2,:).*sol.x);
end
