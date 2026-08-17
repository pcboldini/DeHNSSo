function MLA1 = apply_resonator_wall(MLA1, hr, BF, StabGrid, is_compr)
%APPLY_RESONATOR_WALL  Helmholtz-resonator reactive wall boundary condition (BC) over a mouth.
%
%   MLA1 = apply_resonator_wall(MLA1, hr, BF, StabGrid, is_compr)
%
%   Over the resonator mouth, replaces the wall row of the assembled
%   left-hand-side (LHS) operator MLA1 that normally enforces v' = 0 with the
%   reactive wall relation
%
%       v'|wall  -  Y * p'|wall  =  0,
%
%   which couples the mouth wall-normal velocity to the wall pressure through
%   the resonator admittance Y (the closed-loop resonator response). In the
%   compressible case p' is not a state variable; it is rebuilt from the
%   equation of state, so the relation is imposed on the rho' and T' columns
%   instead. An optional edge taper (hr.taper) smooths v'(x) near the mouth
%   rims without changing the total admittance int Y dx, so it does not add or
%   remove net flux. If hr.uslip is set, the wall u'-row is also replaced by a
%   free-slip condition du'/dy = 0 over the mouth, modeling the fact that the
%   wall is open there.
%
%   Inputs:
%     MLA1      — assembled base LHS operator (modified over the mouth)
%     hr        — resonator settings: Y (admittance), mask (mouth nodes),
%                 optional taper (edge fraction), x (streamwise coordinate),
%                 uslip (free-slip u' over the mouth)
%     BF        — base flow (U, plus rho/T/Ma/gamma for the compressible case)
%     StabGrid  — stability grid (D1 eta-derivative operator, used for u'-slip)
%     is_compr  — true for the (u,v,w,rho,T) path, false for (u,v,w,p)
%   Outputs:
%     MLA1      — operator with the reactive wall (+ optional u'-slip) over the mouth

ny    = size(BF.U, 1);
bsize = (4 + double(is_compr)) * ny;
Y     = hr.Y;
mouth = find(logical(hr.mask(:).'));
tf    = 0;   if isfield(hr,'taper'); tf = hr.taper; end   % edge-taper fraction (0 = sharp mouth edges)
xm    = mouth;                                            % fallback: node index, assumes uniform spacing
if isfield(hr,'x') && ~isempty(hr.x); xm = hr.x(mouth); end
tap   = edge_taper(xm, tf);                               % raised-cosine edge weight, 1 in the interior

% Reactive wall: v'|wall - Y*p'|wall = 0 over the mouth
if is_compr
    ga  = 1.4;  if isfield(BF,'gamma') && ~isempty(BF.gamma); ga = BF.gamma; end
    gM2 = ga * BF.Ma^2;
    for k = 1:numel(mouth)
        i    = mouth(k);   off = bsize*(i-1);
        vrow = off + 2*ny;                  % wall v'-row
        crho = off + 4*ny;                  % wall rho'-column
        cT   = off + 5*ny;                  % wall T'-column
        MLA1(vrow, crho) = -tap(k)*Y * BF.T(ny,i)   / gM2;   % -Y * (dp'/drho') at the wall
        MLA1(vrow, cT)   = -tap(k)*Y * BF.rho(ny,i) / gM2;   % -Y * (dp'/dT')   at the wall
    end
else
    for k = 1:numel(mouth)
        i    = mouth(k);   off = bsize*(i-1);
        vrow = off + 2*ny;                  % wall v'-row
        cp   = off + 4*ny;                  % wall p'-column (p is a state variable)
        MLA1(vrow, cp) = -tap(k)*Y;         % v' - Y*p' = 0
    end
end

% Optional free-slip u' over the mouth (wall opening): replace the u'|wall = 0
% condition with du'/dy = 0. This is applied to MLA1 only, so only the
% fundamental mode sees it.
if isfield(hr,'uslip') && hr.uslip
    D1w = StabGrid.D1(ny, :);               % wall row of the eta-derivative operator
    for k = 1:numel(mouth)
        off  = bsize*(mouth(k)-1);
        urow = off + ny;                    % wall u'-row (u-block comes first, wall is the last eta point)
        MLA1(urow, :)            = 0;        % clear the u' = 0 identity row
        MLA1(urow, off+1:off+ny) = D1w;      % impose du'/dy = 0 (free slip) instead
    end
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
