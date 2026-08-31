function [MLA1, rows] = apply_resonator_wall(MLA1, hr, BF, StabGrid, is_compr)
%APPLY_RESONATOR_WALL  Helmholtz-resonator reactive wall boundary condition (BC) over a mouth.
%
%   [MLA1, rows] = apply_resonator_wall(MLA1, hr, BF, StabGrid, is_compr)
%
%   Over the mouth, replaces the wall row enforcing v' = 0 with
%
%       v'|wall  -  Y * p'|wall  =  0,
%
%   where hr.Y is a scalar (uniform admittance) or one value per mouth node
%   (a profile Y(x)), and hr.mask marks the mouth; setup_boundary resolves both.
%   In the compressible case p' is not a state variable, so
%   the relation is imposed on the rho' and T' columns. hr.taper smooths v'(x)
%   at the rims at fixed int Y dx; hr.uslip frees u' over the mouth.

ny    = size(BF.U, 1);
bsize = (4 + double(is_compr)) * ny;
Y     = hr.Y;
mouth = find(logical(hr.mask(:).'));
tf    = 0;   if isfield(hr,'taper'); tf = hr.taper; end   % edge-taper fraction (0 = sharp mouth edges)
xm    = StabGrid.xun(1, mouth);
if isfield(hr,'x') && ~isempty(hr.x); xm = hr.x(mouth); end
tap   = edge_taper(xm, tf);                               % raised-cosine edge weight, 1 in the interior
if isscalar(Y); Yv = Y * ones(1, numel(mouth)); else; Yv = Y(:).'; end
if numel(Yv) ~= numel(mouth)
    error('apply_resonator_wall:SizeMismatch', ...
          'hr.Y must be a scalar or have %d entries, one per mouth node.', numel(mouth));
end

rows = zeros(1, numel(mouth));

% Reactive wall: v'|wall - Y*p'|wall = 0 over the mouth
if is_compr
    ga  = 1.4;  if isfield(BF,'gamma') && ~isempty(BF.gamma); ga = BF.gamma; end
    gM2 = ga * BF.Ma^2;
    for k = 1:numel(mouth)
        i    = mouth(k);   off = bsize*(i-1);
        vrow = off + 2*ny;                  % wall v'-row
        crho = off + 4*ny;                  % wall rho'-column
        cT   = off + 5*ny;                  % wall T'-column
        MLA1(vrow, crho) = -tap(k)*Yv(k) * BF.T(ny,i)   / gM2;   % -Y * (dp'/drho') at the wall
        MLA1(vrow, cT)   = -tap(k)*Yv(k) * BF.rho(ny,i) / gM2;   % -Y * (dp'/dT')   at the wall
        rows(k) = vrow;
    end
else
    for k = 1:numel(mouth)
        i    = mouth(k);   off = bsize*(i-1);
        vrow = off + 2*ny;                  % wall v'-row
        cp   = off + 4*ny;                  % wall p'-column (p is a state variable)
        MLA1(vrow, cp) = -tap(k)*Yv(k);         % v' - Y*p' = 0
        rows(k) = vrow;
    end
end

% Optional free-slip u' over the mouth (wall opening): replace the u'|wall = 0
% condition with du'/dy = 0. This is applied to MLA1 only, so only the
% fundamental mode sees it.
if isfield(hr,'uslip') && hr.uslip
    D1w  = StabGrid.D1(ny, :);              % wall row of the eta-derivative operator
    urows = zeros(1, numel(mouth));
    for k = 1:numel(mouth)
        off  = bsize*(mouth(k)-1);
        urow = off + ny;                    % wall u'-row (u-block comes first, wall is the last eta point)
        MLA1(urow, :)            = 0;        % clear the u' = 0 identity row
        MLA1(urow, off+1:off+ny) = D1w;      % impose du'/dy = 0 (free slip) instead
        urows(k) = urow;
    end
    rows = [rows, urows];
end
end

function w = edge_taper(xm, frac)
%EDGE_TAPER  Raised-cosine ramp at the mouth edges, on the grid xm.
%   Rescaled so int w dx = int 1 dx over the mouth (via trapezoidal
%   quadrature), so the total admittance stays the same regardless of node
%   spacing.
np = numel(xm);
w  = ones(1, np);
if frac <= 0 || np < 3; return; end
Nt = max(1, round(frac*np));
rr = 0.5 * (1 - cos(pi*(0:Nt-1)/Nt));
w(1:Nt)         = rr;
w(end-Nt+1:end) = fliplr(rr);
xm = xm(:).';
w  = w * (xm(end) - xm(1)) / trapz(xm, w);
end
