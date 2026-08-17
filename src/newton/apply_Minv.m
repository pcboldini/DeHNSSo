function dq_real = apply_Minv(r_real, ctx)
%APPLY_MINV  Apply the preconditioner, working in real arithmetic.
%
%   dq_real = apply_Minv(r_real, ctx)
%
%   Reads r_real as [real(r_phys); imag(r_phys)], applies a complex
%   preconditioner, and returns [real(dq_phys); imag(dq_phys)].
%
%   Three preconditioner modes, chosen by ctx flags:
%     ctx.use_full_precond = true  -> use ctx.LU_full, a single sparse
%       LU/ILU factor of the full coupled Jacobian (see
%       factor_precond_full.m). Uses the most memory but converges best.
%     ctx.use_poly_precond = true  -> polynomial (Neumann-series)
%       preconditioner: per-mode block LU as D^{-1}, ctx.Jh_poly as N.
%       Needs no coupled factorization, so memory grows only linearly
%       with mesh refinement. Number of terms set by ctx.poly_K.
%     otherwise                     -> per-mode block-diagonal solve via
%       ctx.LU{k}. Ignores all coupling between modes: cheapest per call,
%       but convergence stalls at high inflow amplitude A0.

nphys = ctx.nphys;
bsize = ctx.bsize;
N     = nphys * bsize;

% Real-form coupled ILU factor: acts directly on r_real, no complex
% decomposition needed. Captures both the holomorphic and antiholomorphic
% Jacobian parts; used to diagnose the Armijo stall seen with use_full_precond.
if isfield(ctx, 'use_full_real_precond') && ctx.use_full_real_precond
    s = ctx.LU_full_real;
    dq_real = s.U \ (s.L \ r_real);
    return;
end

r_complex = r_real(1:N) + 1i * r_real(N+1:2*N);

if isfield(ctx, 'use_full_precond') && ctx.use_full_precond
    % (I)LU factor of the coupled Jacobian's holomorphic part A only
    % (the antiholomorphic part Ja is ignored here).
    s = ctx.LU_full;
    dq_complex = zeros(N, 1);
    dq_complex(s.q) = s.U \ (s.L \ r_complex(s.p));

elseif isfield(ctx, 'use_poly_precond') && ctx.use_poly_precond
    % Polynomial (Neumann-series) preconditioner:
    %   J = D - Jh, with D = blkdiag(ML), Jh = the cross-mode part of the
    %   NLT (nonlinear-terms) Jacobian
    %   J^{-1} v ≈ sum_{k=0}^{K-1} (D^{-1} Jh)^k D^{-1} v
    K  = ctx.poly_K;
    Jh = ctx.Jh_poly;
    % Term k=0
    x = apply_blockLU(r_complex, ctx);
    dq_complex = x;
    % Terms k=1..K-1
    for kk = 2:K
        Jh_x = Jh * x;
        x    = apply_blockLU(Jh_x, ctx);
        dq_complex = dq_complex + x;
    end

else
    % Block-Jacobi: per-mode only, no coupling between modes.
    dq_complex = apply_blockLU(r_complex, ctx);
end

dq_real = [real(dq_complex); imag(dq_complex)];
end


function dq_complex = apply_blockLU(r_complex, ctx)
%APPLY_BLOCKLU  Per-mode block-diagonal solve using the cached factors.
%   Cached factor type (per mode) is one of:
%     L/U/p/q populated   -> sparse LU / ILU (MATLAB lu/ilu)
%     decomp non-empty    -> LAPACK banded LU (decomposition(A,'banded'))
nphys = ctx.nphys;
bsize = ctx.bsize;
N     = nphys * bsize;
dq_complex = zeros(N, 1);
for k = 1:nphys
    off  = (k-1)*bsize;
    r_k  = r_complex(off+1:off+bsize);
    lu_s = ctx.LU{k};
    if isfield(lu_s, 'decomp') && ~isempty(lu_s.decomp)
        dq_complex(off+1:off+bsize) = lu_s.decomp \ r_k;
    else
        x_k  = zeros(size(r_k));
        x_k(lu_s.q) = lu_s.U \ (lu_s.L \ r_k(lu_s.p));
        dq_complex(off+1:off+bsize) = x_k;
    end
end
end
