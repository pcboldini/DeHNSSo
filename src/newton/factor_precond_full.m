function ctx = factor_precond_full(ctx, q_real)
%FACTOR_PRECOND_FULL  Build and LU-factor the full mode-coupled
%   preconditioner.
%
%   ctx = factor_precond_full(ctx, q_real)
%
%   Assembles M = blkdiag(ML_j) - Jh_NLT(q), where Jh_NLT is the holomorphic
%   part of the analytical Jacobian of the NLT (nonlinear terms), computed
%   by build_JNL. The antiholomorphic part Ja_NLT is dropped from the
%   preconditioner: GMRES still sees the full operator through the
%   Jacobian-free Newton-Krylov (JFNK) matrix-vector product (matvec), so
%   the preconditioner only needs to approximate it. Dropping Ja_NLT keeps
%   memory tractable, since
%   this way the preconditioner is one N x N complex sparse matrix instead
%   of the 2N x 2N real form.
%
%   Sets ctx.LU_full and ctx.use_full_precond = true. After this call,
%   apply_Minv uses this global factorization instead of the per-mode
%   blocks.
%
%   Cost: build_JNL ~5-10 s, assembly ~1 s, sparse LU ~30-90 s.

%% Holomorphic NLT Jacobian
N = ctx.nphys * ctx.bsize;
q_complex = q_real(1:N) + 1i * q_real(N+1:2*N);

% Free the previous global factor before rebuilding.
if isfield(ctx, 'LU_full'); ctx.LU_full = []; end

t_jnl = tic;
[Jh_NLT, ~] = build_JNL(q_complex, ctx);
fprintf('   factor_full: build_JNL %.1fs (nnz=%.1fM)\n', toc(t_jnl), nnz(Jh_NLT)/1e6);

% Drop tiny entries (below 1e-14 of the max magnitude). This shaves ~20-30%
% off the nonzero count without hurting the preconditioner's quality.
mx = max(abs(nonzeros(Jh_NLT)));
if mx > 0
    Jh_NLT = Jh_NLT .* (abs(Jh_NLT) > 1e-14*mx);
end

%% Assemble M = blkdiag(ML_j) - Jh_NLT
% Free ML_block as soon as the sum is formed.
t_asm = tic;
ML_block = blkdiag(ctx.ML{:});
M_full = ML_block - Jh_NLT;
clear ML_block Jh_NLT
fprintf('   factor_full: assemble %.1fs (nnz(M)=%.1fM)\n', toc(t_asm), nnz(M_full)/1e6);

% Free the per-mode LU cache to recover several GB for the global factor.
ctx.LU = cell(ctx.nphys, 1);

%% Incomplete-LU factorization
% Use incomplete LU (ILU) as the preconditioner: memory scales with nnz(M)
% and it builds in seconds. Try ILU(0) first (no fill-in); if that breaks
% down, fall back to a thresholded, pivoted variant, which is more robust
% at the cost of some extra fill-in.
t_lu = tic;
L = []; U = [];
try
    setup0 = struct('type','nofill');
    [L, U] = ilu(M_full, setup0);
    kind = 'ilu0';
catch ME0
    fprintf('   factor_full: ilu(0) failed (%s), trying ilutp ...\n', ME0.message);
end
if isempty(L)
    try
        setup1 = struct('type','ilutp','droptol',1e-4,'thresh',0.1);
        [L, U] = ilu(M_full, setup1);
        kind = 'ilutp';
    catch ME1
        error('factor_precond_full:iluFailed', ...
              'Both ILU(0) and ILUTP(1e-4) broke down: %s', ME1.message);
    end
end
fprintf('   factor_full: %s %.1fs (L+U: %.0f MB)\n', kind, toc(t_lu), (nnz(L)+nnz(U))*16/1e6);
clear M_full

%% Store factor for apply_Minv
% Identity row/column permutations, so apply_Minv can use the same generic
% dispatch logic (it applies p on the input and q on the output).
N_total = ctx.nphys * ctx.bsize;
id = (1:N_total).';
ctx.LU_full = struct('L', L, 'U', U, 'p', id, 'q', id);
ctx.use_full_precond = true;
end
