function BF_eff = build_effective_BF(BF_orig, q_complex, ctx)
%BUILD_EFFECTIVE_BF  Augment base flow with the current MFD perturbation.
%
%   BF_eff = build_effective_BF(BF_orig, q_complex, ctx)
%
%   Returns BF_eff with U,V,W and dxU,dyU,... (etc.) augmented by the
%   real part of the current mean-flow-distortion (MFD) state. Calling
%   LHS(BF_eff, ...) then produces per-mode linear operators that include
%   the linearisation of the dominant nonlinear self-coupling: each
%   physical mode being convected by the MFD, plus the local "stretching"
%   terms from the MFD's gradients.
%
%   This is the DeHNSSo equivalent of Scholten et al. (2024)'s "simplified
%   Jacobian without cross-derivative terms": the cross-mode (higher
%   harmonic) couplings are dropped (still iterated by Newton), but the
%   block-diagonal MFD self-coupling is retained — that's where the
%   stiffness lives at large amplitudes.
%
%   By symmetry, MFD is real (StabRes.q(j_MFD,:) is its own conjugate),
%   so we take real() defensively to handle imaginary roundoff during
%   Newton iterations.

phys  = ctx.phys_modes;
nphys = ctx.nphys;
ny    = ctx.ny;
nx    = ctx.nx;
bsize = ctx.bsize;

% Locate MFD mode (omega=0, beta=0) within physical modes
k_MFD = [];
for k = 1:nphys
    j = phys(k);
    if ctx.omegavec(j) == 0 && ctx.betavec(j) == 0
        k_MFD = k; break
    end
end

BF_eff = BF_orig;
if isempty(k_MFD); return; end

% Extract MFD perturbation
off  = (k_MFD-1)*bsize;
qmfd = reshape(q_complex(off+1:off+bsize), 4*ny, nx);
u_mfd = real(qmfd(1:ny,        :));
v_mfd = real(qmfd(ny+1:2*ny,   :));
w_mfd = real(qmfd(2*ny+1:3*ny, :));

% Physical derivatives via metric chain rule.
%   ∂/∂x = ξx·∂/∂ξ + ηx·∂/∂η   (and similarly for ∂/∂y)
% ξ-derivative uses the sparse non-uniform FD matrix Grid.D_xi (standard
% sign convention; Fornberg-based). Reduces to 4th-order uniform FD on
% uniform grids, bit-identical to the old -FD1d4o(·, dxi) path there.
D_xi   = ctx.D_xi;
D1     = ctx.D1;
xix    = ctx.xix;   etax = ctx.etax;
xiy    = ctx.xiy;   etay = ctx.etay;

dudxi  = u_mfd * D_xi.';
dvdxi  = v_mfd * D_xi.';
dwdxi  = w_mfd * D_xi.';
dudeta = D1 * u_mfd;
dvdeta = D1 * v_mfd;
dwdeta = D1 * w_mfd;

BF_eff.U   = BF_orig.U   + u_mfd;
BF_eff.V   = BF_orig.V   + v_mfd;
BF_eff.W   = BF_orig.W   + w_mfd;
BF_eff.dxU = BF_orig.dxU + dudxi.*xix + dudeta.*etax;
BF_eff.dxV = BF_orig.dxV + dvdxi.*xix + dvdeta.*etax;
BF_eff.dxW = BF_orig.dxW + dwdxi.*xix + dwdeta.*etax;
BF_eff.dyU = BF_orig.dyU + dudxi.*xiy + dudeta.*etay;
BF_eff.dyV = BF_orig.dyV + dvdxi.*xiy + dvdeta.*etay;
BF_eff.dyW = BF_orig.dyW + dwdxi.*xiy + dwdeta.*etay;

end
