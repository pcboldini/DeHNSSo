function [dxx,dyy,dxy] = phys_d2(q, qxi, qeta, D_xi2, D1, D2, m)
%PHYS_D2  Curvilinear second physical derivatives of a field q(ny,nx).
%
%   [dxx, dyy, dxy] = phys_d2(q, qxi, qeta, D_xi2, D1, D2, m)
%
%   Applies the full curvilinear chain rule (first AND second inverse metrics),
%   matching gridgen/compute_bf_derivatives:
%     d_xx q = xix^2 q_xixi + 2 xix etax q_xieta + etax^2 q_etaeta + xixx q_xi + etaxx q_eta
%     d_yy q = xiy^2 q_xixi + 2 xiy etay q_xieta + etay^2 q_etaeta + xiyy q_xi + etayy q_eta
%     d_xy q = xix xiy q_xixi + (xix etay + xiy etax) q_xieta + etax etay q_etaeta + xixy q_xi + etaxy q_eta
%
%   Inputs:
%     q       — field on the (ny x nx) grid
%     qxi     — computational first derivative q*D_xi.' (streamwise)
%     qeta    — computational first derivative D1*q (wall-normal)
%     D_xi2   — streamwise second-derivative operator (4th-order)
%     D1, D2  — wall-normal first / second derivative operators
%     m       — inverse metrics (xix,etax,xiy,etay,xixx,etaxx,xiyy,etayy,xixy,etaxy)
%   Outputs:
%     dxx, dyy, dxy — second physical derivatives (dxy only if requested)

% Computational second derivatives
qxixi   = q * D_xi2.';        % d^2/dxi^2  (4th-order streamwise operator)
qetaeta = D2 * q;            % d^2/deta^2 (wall-normal spectral/FD operator)
qxieta  = D1 * qxi;         % d/deta of d/dxi

% Assemble physical derivatives via the curvilinear chain rule
dxx = m.xix.^2.*qxixi + 2*m.xix.*m.etax.*qxieta + m.etax.^2.*qetaeta + m.xixx.*qxi + m.etaxx.*qeta;
dyy = m.xiy.^2.*qxixi + 2*m.xiy.*m.etay.*qxieta + m.etay.^2.*qetaeta + m.xiyy.*qxi + m.etayy.*qeta;
if nargout > 2
    dxy = m.xix.*m.xiy.*qxixi + (m.xix.*m.etay + m.xiy.*m.etax).*qxieta ...
        + m.etax.*m.etay.*qetaeta + m.xixy.*qxi + m.etaxy.*qeta;
end
end
