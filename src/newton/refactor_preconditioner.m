function ctx = refactor_preconditioner(ctx, q_real)
%REFACTOR_PRECONDITIONER  Rebuild the LU preconditioner using current MFD.
%
%   Updates ctx.LU (the per-mode preconditioner factorisation) by calling
%   LHS.m with an effective base flow that includes the current
%   mean-flow-distortion perturbation. Original ctx.ML matrices are NOT
%   touched — they are used by residual_hns to evaluate F(q) and must
%   remain the true linear operator.
%
%   Cost: one full LHS rebuild (~2-3 s) + nphys LU factorisations
%   (~4 s each, so ~20-25 s for a 5-mode 2D run). Free old factors first
%   to keep peak memory bounded.

N = ctx.nphys * ctx.bsize;
q_complex = q_real(1:N) + 1i * q_real(N+1:2*N);

BF_eff = build_effective_BF(ctx.BF_orig, q_complex, ctx);

t_lhs = tic;
[MLAMFD_e, MLA1_e, MLA2_e, MLA3_e, MLA4_e] = ...
    LHS(BF_eff, ctx.Bound, ctx.StabGrid, ctx.Init, ctx.Stab, ctx.StabRes, ctx.Opt);
MFDdiff_e = MLAMFD_e - MLA1_e;
clear MLAMFD_e
fprintf('   refactor: LHS rebuilt in %.1f s\n', toc(t_lhs));

% Free old LU factors before allocating new ones (peak memory control)
ctx.LU = cell(ctx.nphys, 1);

key_cache = struct();
for k = 1:ctx.nphys
    j     = ctx.phys_modes(k);
    omega = ctx.omegavec(j);
    beta  = ctx.betavec(j);
    key   = ml_key_local(omega, beta);

    if ~isfield(key_cache, key)
        if omega == 0 && beta == 0
            ML = MLA1_e + MFDdiff_e;
        elseif omega == 0
            ML = MLA1_e + MLA3_e*beta + MLA4_e*beta^2;
        elseif beta == 0
            ML = MLA1_e + MLA2_e*omega;
        else
            ML = MLA1_e + MLA2_e*omega + MLA3_e*beta + MLA4_e*beta^2;
        end
        tic_lu = tic;
        [L, U, p, q_p] = lu(ML, 'vector');
        fprintf('   refactor: LU (%d,%d) %.1fs\n', ...
                ctx.Mvec(j), ctx.Nvec(j), toc(tic_lu));
        key_cache.(key) = struct('L', L, 'U', U, 'p', p, 'q', q_p);
        clear ML L U
    end
    ctx.LU{k} = key_cache.(key);
end
end

function key = ml_key_local(omega, beta)
if omega == 0 && beta == 0
    key = 'MFD';
else
    key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
end
end
