function ctx = build_ctx(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes)
%BUILD_CTX  Assemble Newton solver context: LHS, per-mode LU, BC vectors.
%
%   ctx = build_ctx(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes)
%
%   Assembles the per-mode linear operators (LHS), their block preconditioner
%   factors, and the inhomogeneous BC vectors into a single context struct
%   consumed by the Newton-Krylov / pseudo-arclength solver.
%
%   Used by solve_newton.m (pseudo-arclength continuation in A0_fund) and by
%   the build_JNL instrumentation in paper/checks/, which call it directly and
%   set newton_precond='none' to skip the per-mode block LUs they don't need.

%% Assemble LHS operators
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
nvar     = 4 + is_compr;   % 4 incompressible | 5 compressible
fprintf(' building LHS (%s) ...', ternary_c(is_compr,'compressible','incompressible'));
t = tic;
MLA5 = [];
if is_compr
    [MLAMFD, MLA1, MLA2, MLA3, MLA4, MLA5] = LHS_compressible(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
else
    [MLAMFD, MLA1, MLA2, MLA3, MLA4]       = LHS(BF, Bound, StabGrid, Init, Stab, StabRes, Opt);
end
MFD_diff = MLAMFD - MLA1;
clear MLAMFD

% Cheap-stencil LHS matrices for the preconditioner (Liu/Brehm pattern).
% When Opt.newton_block_factor='cheap', the per-mode preconditioner is the
% full LU of a 2nd-order-FD Jacobian (cheap), while the Newton matvec keeps
% using the full-order operator from MLA1..MLA4 above.
use_cheap_precond = isfield(Opt,'newton_block_factor') && ...
                    strcmpi(Opt.newton_block_factor, 'cheap') && ~is_compr;
if is_compr && isfield(Opt,'newton_block_factor') && strcmpi(Opt.newton_block_factor,'cheap')
    fprintf(' [compressible: cheap-stencil precond not yet supported, using full-order LU] ');
end
if use_cheap_precond
    fprintf(' building cheap-stencil LHS for preconditioner ...');
    t_cheap_lhs = tic;
    StabGrid_cheap = make_cheap_stabgrid(StabGrid);
    [MLAMFD_c, MLA1_c, MLA2_c, MLA3_c, MLA4_c] = ...
        LHS(BF, Bound, StabGrid_cheap, Init, Stab, StabRes, Opt);
    MFD_diff_c = MLAMFD_c - MLA1_c;
    clear MLAMFD_c
    fprintf(' %.1fs\n', toc(t_cheap_lhs));
end
fprintf(' %.1fs\n', toc(t));

%% Active modes and context fields
half  = round(Mode.nf/2);
phys  = Mode.RunJ(Mode.RunJ >= half);   % active physical modes only
nphys = numel(phys);

ctx.phys_modes = phys;
ctx.nphys      = nphys;
ctx.nf         = Mode.nf;
ctx.ny         = StabGrid.ny;
ctx.nx         = StabGrid.nx;
ctx.nvar       = nvar;
ctx.bsize      = nvar*ctx.ny*ctx.nx;
ctx.StabGrid   = StabGrid;
ctx.D1         = StabGrid.D1;
ctx.HB         = Stab.HB;
ctx.omegavec   = StabRes.omegavec;
ctx.betavec    = StabRes.betavec;
ctx.Bound      = Bound;
ctx.Mode       = Mode;
ctx.BF_orig    = BF;
ctx.Init       = Init;
ctx.Stab       = Stab;
ctx.StabRes    = StabRes;
ctx.Opt        = Opt;
ctx.Mvec       = Stab.Mvec;
ctx.Nvec       = Stab.Nvec;
ctx.dxi        = StabGrid.x(2) - StabGrid.x(1);   % legacy scalar (for display only)
ctx.D_xi       = StabGrid.D_xi;                    % sparse non-uniform FD matrix
ctx.xix        = StabGrid.xix;
ctx.etax       = StabGrid.etax;
ctx.xiy        = StabGrid.xiy;
ctx.etay       = StabGrid.etay;

ctx.ML  = cell(nphys, 1);
ctx.LU  = cell(nphys, 1);
ctx.rBC = cell(nphys, 1);

Opt_rhs = Opt;
Opt_rhs.AF = 1;   % store UNIT amplitude BC; residual_hns scales by ctx.AF
key_cache = struct();

% 'none' = don't build the per-mode block LUs at all (only ctx.ML is needed).
% Used by the build_JNL instrumentation in paper/checks/.
skip_block_lu = isfield(Opt,'newton_precond') ...
                && ismember(lower(Opt.newton_precond), {'full','full_real','none'});

% Shared block preconditioner: build one exact banded factor for the wave modes
% + one for the MFD mode and point all modes at them, instead of one factor per
% mode. The modes differ only by omega*MLA2 (a small perturbation), so a single
% exact factor preconditions each block well while cutting the factor memory
% from nphys x to ~2x. ctx.ML{k} stay per-mode (correct operators for the matvec).
share_pc = isfield(Opt,'newton_share_precond') && strcmpi(Opt.newton_share_precond,'on') ...
           && ~skip_block_lu;

%% Per-mode LHS blocks and block preconditioners
for k = 1:nphys
    j     = phys(k);
    omega = StabRes.omegavec(j);
    beta  = StabRes.betavec(j);
    key   = ml_key(omega, beta);

    if ~isfield(key_cache, key)
        % MLA5 = i*beta*d_x cross-viscous block (compressible only; empty otherwise).
        b5 = 0; if is_compr && ~isempty(MLA5); b5 = beta; end
        if omega == 0 && beta == 0
            ML = MLA1 + MFD_diff;
        elseif omega == 0
            ML = MLA1 + MLA3*beta + MLA4*beta^2;
            if b5, ML = ML + MLA5*b5; end
        elseif beta == 0
            ML = MLA1 + MLA2*omega;
        else
            ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2;
            if b5, ML = ML + MLA5*b5; end
        end
        if skip_block_lu || share_pc
            if share_pc
                fprintf('   ML mode (%d,%d) built (shared block-LU precond)\n', ...
                        Stab.Mvec(j), Stab.Nvec(j));
            else
                fprintf('   ML mode (%d,%d) built (no block-LU, full ILU will be used)\n', ...
                        Stab.Mvec(j), Stab.Nvec(j));
            end
            key_cache.(key) = struct('ML', ML, 'L', [], 'U', [], 'p', [], 'q', [], 'decomp', []);
        else
            block_factor = 'lu';
            if isfield(Opt,'newton_block_factor') && ~isempty(Opt.newton_block_factor)
                block_factor = lower(Opt.newton_block_factor);
            end
            decomp = [];
            if strcmp(block_factor, 'ilu0')
                fprintf('   ILU0 mode (%d,%d) ...', Stab.Mvec(j), Stab.Nvec(j));
                t_lu = tic;
                [L, U] = ilu(ML, struct('type','nofill'));
                idn = (1:size(ML,1)).';
                p   = idn;
                q_p = idn;
                fprintf(' %.1fs (%.0f MB)\n', toc(t_lu), (nnz(L)+nnz(U))*16/1e6);
            elseif strcmp(block_factor, 'banded')
                fprintf('   banded-LU mode (%d,%d) ...', Stab.Mvec(j), Stab.Nvec(j));
                t_lu   = tic;
                decomp = decomposition(ML, 'banded');
                t_dec  = toc(t_lu);
                [il, jl] = find(ML);
                kl_est = max(il - jl);
                ku_est = max(jl - il);
                band_mb = (2*kl_est + ku_est + 1) * size(ML,1) * 16 / 1e6;
                fprintf(' %.1fs (band %.0f MB, KL=%d KU=%d)\n', ...
                        t_dec, band_mb, kl_est, ku_est);
                L = []; U = [];
                idn = (1:size(ML,1)).';
                p   = idn;
                q_p = idn;
            elseif strcmp(block_factor, 'cheap')
                if omega == 0 && beta == 0
                    ML_pre = MLA1_c + MFD_diff_c;
                elseif omega == 0
                    ML_pre = MLA1_c + MLA3_c*beta + MLA4_c*beta^2;
                elseif beta == 0
                    ML_pre = MLA1_c + MLA2_c*omega;
                else
                    ML_pre = MLA1_c + MLA2_c*omega + MLA3_c*beta + MLA4_c*beta^2;
                end
                fprintf('   cheap-LU mode (%d,%d) ...', Stab.Mvec(j), Stab.Nvec(j));
                t_lu = tic;
                [L, U, p, q_p] = lu(ML_pre, 'vector');
                fprintf(' %.1fs (%.0f MB, full-ML nnz=%.1fM cheap-ML nnz=%.1fM)\n', ...
                    toc(t_lu), (nnz(L)+nnz(U))*16/1e6, nnz(ML)/1e6, nnz(ML_pre)/1e6);
            else
                fprintf('   LU mode (%d,%d) ...', Stab.Mvec(j), Stab.Nvec(j));
                t_lu = tic;
                [L, U, p, q_p] = lu(ML, 'vector');
                fprintf(' %.1fs (%.0f MB)\n', toc(t_lu), (nnz(L)+nnz(U))*16/1e6);
            end
            key_cache.(key) = struct('ML', ML, 'L', L, 'U', U, ...
                                     'p', p, 'q', q_p, 'decomp', decomp);
        end
    end
    entry = key_cache.(key);
    ctx.ML{k} = entry.ML;
    if ~skip_block_lu
        ctx.LU{k} = struct('L', entry.L, 'U', entry.U, ...
                           'p', entry.p, 'q', entry.q, 'decomp', entry.decomp);
    end
    ctx.rBC{k} = RHS(Bound, StabGrid, Mode, Opt_rhs, StabRes, j);
end

%% Shared block preconditioner (optional)
if share_pc
    % Build (at most) two shared banded factors: one for the MFD mode, one for
    % the wave modes. Point every ctx.LU{k} at the matching factor (MATLAB
    % copy-on-write keeps them sharing memory since apply_blockLU never mutates).
    fprintf('   building shared banded factor(s) ...');
    t_sh = tic;  mfd_fac = [];  wave_fac = [];
    for k = 1:nphys
        j = phys(k);  om = StabRes.omegavec(j);  bt = StabRes.betavec(j);
        if om==0 && bt==0
            if isempty(mfd_fac);  mfd_fac  = decomposition(ctx.ML{k}, 'banded'); end
        else
            if isempty(wave_fac); wave_fac = decomposition(ctx.ML{k}, 'banded'); end
        end
    end
    if isempty(wave_fac); wave_fac = mfd_fac; end
    if isempty(mfd_fac);  mfd_fac  = wave_fac; end
    for k = 1:nphys
        j = phys(k);  om = StabRes.omegavec(j);  bt = StabRes.betavec(j);
        if om==0 && bt==0, fac = mfd_fac; else, fac = wave_fac; end
        ctx.LU{k} = struct('L', [], 'U', [], 'p', [], 'q', [], 'decomp', fac);
    end
    fprintf(' %.1fs (<=2 factors shared across %d modes)\n', toc(t_sh), nphys);
end
end


function key = ml_key(omega, beta)
if omega == 0 && beta == 0
    key = 'MFD';
else
    key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
end
end

function s = ternary_c(c,a,b); if c, s=a; else, s=b; end; end
