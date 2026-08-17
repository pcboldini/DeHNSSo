function [StabRes, info] = solve_newton(Stab, StabGrid, Opt)
%SOLVE_NEWTON  Pseudo-arclength continuation in the inflow amplitude A0_fund.
%
%   [StabRes, info] = solve_newton(Stab, StabGrid, Opt)
%
%   Continues the nonlinear HNS solution along the branch A0_fund = lambda,
%   from Opt.newton_start up to Opt.newton_target, and returns the converged state at
%   the target. Use when a direct solve at the target amplitude will not
%   converge from a cold start (high-amplitude K-type / Scholten cases).
%
%   Method (Keller pseudo-arclength + Shamanskii corrector):
%     predictor  — tangent + secant curvature extrapolation in (q, lambda)
%     corrector  — modified Newton on the full coupled real-form Jacobian,
%                  factorised once with MUMPS and REUSED across steps
%                  (back-solves only) until it ages out or stalls.
%
%   The reuse is what makes this fast: the converged-state factorisation costs
%   ~230x a back-solve, and along the branch the state moves little between
%   steps, so one factor stays valid for several steps. Refactor triggers:
%   age >= Opt.newton_reuse, or the previous step needed > Opt.newton_it_thresh
%   corrector iterations, or the corrector fails (stale-factor fallback).
%
%   lambda enters F only through the inflow Dirichlet BC, hence affinely:
%       F(q,lambda) = residual_hns(q,ctx0) + (lambda-lam0)*F_lambda
%   F_lambda is obtained once by finite difference and its linearity is
%   checked and reported.
%
%   Opt fields (defaults):
%     newton_target     — target A0_fund (required)
%     newton_start       — continuation start   (default Stab.A0_fund)
%     newton_seed2       — second seed point    (default 1.04*newton_start)
%     newton_reuse       — reuse steps before forced refactor (default 8)
%     newton_it_thresh  — refactor next step if corrector exceeded this (default 8)
%     newton_maxit      — max corrector back-solves before stale (default 25)
%     newton_dsmax      — max arclength step (default 3)
%     newton_maxstep    — max continuation steps (default 400)
%     newton_tol        — corrector residual tolerance (default 1e-6)
%     newton_seedfile   — cache for the two seed points ([] = no cache)
%     newton_ckfile     — checkpoint file written each step ([] = no checkpoint)
%     newton_resume     — 'on' to resume from newton_ckfile (default 'off')
%
%   See also BUILD_CTX, BUILD_JNL, RESIDUAL_HNS.

%% Options and continuation parameters
t_all = tic; cum = 0;

if ~isfield(Opt,'newton_target') || isempty(Opt.newton_target)
    error('solve_newton:noTarget','Opt.newton_target (target A0_fund) is required.');
end
A_target = Opt.newton_target;
lam0 = Stab.A0_fund; if isfield(Opt,'newton_start') && ~isempty(Opt.newton_start); lam0 = Opt.newton_start; end
lam1 = 1.04*lam0;    if isfield(Opt,'newton_seed2') && ~isempty(Opt.newton_seed2); lam1 = Opt.newton_seed2; end
Kmax        = get_o(Opt,'newton_reuse', 8);
it_thresh   = get_o(Opt,'newton_it_thresh', 8);
reuse_maxit = get_o(Opt,'newton_maxit', 25);
dsmax       = get_o(Opt,'newton_dsmax', 3);
maxstep     = get_o(Opt,'newton_maxstep', 400);
arctol      = get_o(Opt,'newton_tol', 1e-6);
seedfile    = get_o(Opt,'newton_seedfile', []);
ckfile      = get_o(Opt,'newton_ckfile', []);
do_resume   = isfield(Opt,'newton_resume') && strcmpi(Opt.newton_resume,'on');

fprintf('\n--- Pseudo-arclength continuation: A0 %.4e -> %.4e (u_rms %.4f) ---\n', ...
        lam0, A_target, sqrt(2)*A_target);

% Base flow fields for build_ctx (as assembled by the callers/drivers).
BF.Re  = StabGrid.Re;
BF.U   = StabGrid.U;   BF.V   = StabGrid.V;   BF.W   = StabGrid.W;
BF.dxU = StabGrid.dxU; BF.dxV = StabGrid.dxV; BF.dxW = StabGrid.dxW;
BF.dyU = StabGrid.dyU; BF.dyV = StabGrid.dyV; BF.dyW = StabGrid.dyW;

% Compressible base flow requires the extra thermodynamic/transport fields
% (LHS_compressible reads BF.Ma, rho, T and their derivatives). Mirror main.m.
if isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on')
    req = {'rho','T','dxrho','dyrho','dxT','dyT','dxxU','dxxV','dxxW','dxxT', ...
           'ddyU','ddyV','ddyW','ddyT','dxyU','dxyV', ...
           'Ma','Pr','gamma','Sstar','mu','dmuT','d2muT2'};
    for kf = 1:numel(req); BF.(req{kf}) = StabGrid.(req{kf}); end
end

%% Context at a given lambda
% Names here are deliberately disjoint from the parent workspace: this is a
% nested function, so any name shared with the parent is the SAME variable
% (an internal `qc`/`q0`/`Mode` would clobber the continuation state).
    function [ctx_o, qv_o, SG_o, Bnd_o, Mode_o, Opt_o, Stab_o, SR_o] = mkctx(lam)
        S_l = Stab; S_l.A0_fund = lam;
        [SG_o, Init_l, Mode_o, Opt_o, Stab_o, SR_o] = init(StabGrid, Opt, S_l);
        [Mode_o, SR_o] = inflow(Mode_o, Stab_o, SG_o, SR_o, Opt_o);
        [Bnd_o]        = setup_boundary(SG_o, Opt_o, Stab_o);
        Mode_o.RunJ    = 1:Mode_o.nf; Mode_o.Aactive = ones(1, Mode_o.nf);
        ctx_o = build_ctx(BF, Bnd_o, SG_o, Init_l, Mode_o, Opt_o, Stab_o, SR_o);
        ctx_o.AF = 1;
        Nl_l = ctx_o.nphys*ctx_o.bsize; qc_l = zeros(Nl_l,1);
        for k_l = 1:ctx_o.nphys
            j_l = ctx_o.phys_modes(k_l); off_l = (k_l-1)*ctx_o.bsize;
            qc_l(off_l+1:off_l+ctx_o.bsize) = SR_o.q(j_l,:).';
        end
        qv_o = [real(qc_l); imag(qc_l)];
    end

fprintf(' building ctx0 @ %.4e ...\n', lam0);
[ctx0, q_init, SG0, Bnd0, Mode0, Opt0, Stab0, SR0] = mkctx(lam0);
N = ctx0.nphys*ctx0.bsize; twoN = 2*N;
ML = blkdiag(ctx0.ML{:});
fprintf(' N(complex)=%d  2N=%d  nphys=%d bsize=%d nf=%d\n', ...
        N, twoN, ctx0.nphys, ctx0.bsize, ctx0.nf);

% MUMPS ordering must follow the wall-normal discretisation, or it fails with
% INFOG(1)=-13 on the other operator (measured Jul 2026, lsl, filled M=4):
%   Chebyshev (dense D1, ny-cliques)   AMD  123 GB ok / METIS 216 GB -13
%   FD        (sparse D1, 2D stencil)  METIS 181 GB ok / AMD   311 GB -13
% Detect by whether the wall-normal operator is sparse (FD) or dense (Cheb),
% and set the default that mumps_factor reads via get_o. Explicit
% Opt.newton_mumps_ordering still overrides.
if ~isfield(Opt,'newton_mumps_ordering') || isempty(Opt.newton_mumps_ordering)
    is_fd = isfield(ctx0,'D1') && issparse(ctx0.D1);
    Opt.newton_mumps_ordering = tern(is_fd, 7, 0);   % 7=auto(->METIS on FD), 0=AMD
    fprintf(' MUMPS ordering: %s (ICNTL(7)=%d), auto-selected from %s wall-normal operator\n', ...
            tern(is_fd,'auto/METIS','AMD'), Opt.newton_mumps_ordering, tern(is_fd,'sparse FD','dense Cheb'));
end

%% dF/dlambda (affine in lambda: inflow BC only)
h = 1e-6; [ctxh, ~] = mkctx(lam0 + h);
Fa = residual_hns(q_init, ctx0); Fb = residual_hns(q_init, ctxh);
Flam = (Fb - Fa)/h; phi = -Flam; phi = phi(:); Flam = Flam(:);
[ctx2, ~] = mkctx(lam0 + 50*h); Fchk = residual_hns(q_init, ctx2);
fprintf(' ||F_lambda||=%.3e  linearity rel.err=%.2e\n', ...
        norm(Flam), norm(Fchk - (Fa + 50*h*Flam))/max(norm(Fchk), eps));
clear ctxh ctx2

evalF = @(q,lam) residual_hns(q, ctx0) + (lam - lam0)*Flam;

    function Jr = Jreal(q)
        qcx = q(1:N) + 1i*q(N+1:twoN);
        [Jh, Ja] = build_JNL(qcx, ctx0); A = ML - Jh; B = -Ja;
        Jr = [real(A)+real(B), -imag(A)+imag(B); imag(A)+imag(B), real(A)-real(B)];
    end
% Natural-parameter Newton for the two seed points. All locals are suffixed
% _n: this is a nested function, so an unsuffixed `id`/`Jr`/`F` would BE the
% parent's held MUMPS instance/Jacobian. natnewton destroys its instance each
% iteration, so sharing `id` would leave the parent holding a destroyed handle
% and the first cross-step back-solve would run on freed memory.
    function q_out = natnewton(q_in, lam, tol, mx)
        q_out = q_in;
        for it_n = 1:mx
            F_n = evalF(q_out, lam); nF_n = norm(F_n);
            fprintf('   nat it%d lam=%.5e ||F||=%.3e\n', it_n, lam, nF_n);
            if nF_n < tol, break; end
            Jr_n = Jreal(q_out); id_n = mumps_factor(Jr_n, Opt);
            [dq_n, id_n] = mumps_bsolve(id_n, Jr_n, -F_n);
            id_n = mumps_destroy(id_n); %#ok<NASGU>
            q_out = q_out + dq_n;
        end
    end

%% Two seed points
% Skipped entirely when resuming: the seeds exist only to build the initial
% tangent, and the resume branch below discards q0/q1 and recomputes the
% tangent from the checkpointed state. Running them anyway costs two cold-start
% natural-parameter solves at the resume amplitude (the regime the continuation
% exists to avoid), which on a refined grid is hours of wasted factorisations.
have_ck = do_resume && ~isempty(ckfile) && exist(ckfile,'file') == 2;
if have_ck
    fprintf('\n -- seed points skipped (resuming from %s) --\n', ckfile);
    q0 = q_init; q1 = q_init; tq = zeros(size(q_init)); tl = 1; ds = 0; %#ok<NASGU>
    qc = q_init(:); lc = lam0;
else
    fprintf('\n -- seed points --\n'); t_seed = tic;
    if ~isempty(seedfile) && exist(seedfile,'file') == 2
        L = load(seedfile); q0 = L.q0(:); q1 = L.q1(:);
        fprintf(' loaded seed from %s\n', seedfile);
    else
        q0 = natnewton(q_init, lam0, 1e-6, 30);
        q1 = natnewton(q0,     lam1, 1e-6, 30);
        if ~isempty(seedfile)
            save(seedfile, 'q0', 'q1', 'lam0', 'lam1', '-v7.3');
            fprintf(' saved seed to %s\n', seedfile);
        end
    end
    fprintf(' [timing] seed: %.1f min\n', toc(t_seed)/60);

    dq = q1 - q0; dl = lam1 - lam0; nrm = sqrt(dq'*dq + dl^2);
    tq = dq/nrm; tl = dl/nrm; qc = q1(:); lc = lam1; ds = max(nrm, 3e-5); tq = tq(:);
    fprintf(' initial tangent: dA0/ds=%.3e  ds=%.3e\n', tl, ds);
end

%% Cross-step reuse state
id = []; Jr = []; age = 1e9; it_prev = 0; nfac_total = 0;
traj = []; tq_old = []; tl_old = []; ds_prev = 0.35;

if do_resume && ~isempty(ckfile) && exist(ckfile,'file') == 2
    L = load(ckfile); qc = L.q(:); lc = L.l; traj = L.traj;
    if isfield(L,'nfac_total'); nfac_total = L.nfac_total; end
    Jr = Jreal(qc); id = mumps_factor(Jr, Opt); nfac_total = nfac_total + 1; age = 0;
    [w, id] = mumps_bsolve(id, Jr, phi); tl = 1/sqrt(w'*w + 1); tq = tl*w; tq = tq(:);
    if tl < 0, tl = -tl; tq = -tq; end
    ds = get_o(Opt,'newton_ds', 2.0);
    fprintf(' RESUMED from %s: A0=%.6e u_rms=%.5e, %d steps done, dA0/ds=%.3e, ds=%.3e\n', ...
            ckfile, lc, sqrt(2)*lc, size(traj,1), tl, ds);
end

info.status = 'maxstep'; info.nfac_total = 0; info.steps = 0;

%% Main continuation loop
for step = (size(traj,1)+1):maxstep
    t_step = tic;
    if tl > 0 && lc + ds*tl > A_target, ds = (A_target - lc)/tl; end
    if ~isempty(tq_old)
        cq = (tq - tq_old)/ds_prev; cl = (tl - tl_old)/ds_prev;
        qp = qc + ds*tq + 0.5*ds^2*cq; lp = lc + ds*tl + 0.5*ds^2*cl;
    else
        qp = qc + ds*tq; lp = lc + ds*tl;
    end
    q = qp(:); l = lp; conv = false; nfac_step = 0;

    % reuse the standing factor, or refactor at the predictor
    need_fac = isempty(id) || age >= Kmax || it_prev > it_thresh;
    if need_fac
        if ~isempty(id), id = mumps_destroy(id); id = []; end
        Jr = Jreal(q); id = mumps_factor(Jr, Opt);
        nfac_step = 1; nfac_total = nfac_total + 1; age = 0;
    end
    [b, id] = mumps_bsolve(id, Jr, phi); b = b(:);

    for it = 1:reuse_maxit
        q = q(:); F = evalF(q, l); F = F(:);
        Nc  = tq(:).'*(q - qc) + tl*(l - lc) - ds;
        res = sqrt(F'*F + Nc^2);
        if res < arctol, conv = true; break; end
        [a, id] = mumps_bsolve(id, Jr, -F); a = a(:);
        dl = (-Nc - tq(:).'*a)/(tl + tq(:).'*b); dl = dl(1);
        dq = a + dl*b; q = q(:) + dq(:); l = l + dl;
    end

    % stale-factor fallback: refactor at current q, finish with a fresh factor
    if ~conv
        if ~isempty(id), id = mumps_destroy(id); end
        Jr = Jreal(q); id = mumps_factor(Jr, Opt);
        nfac_step = nfac_step + 1; nfac_total = nfac_total + 1; age = 0;
        [b, id] = mumps_bsolve(id, Jr, phi); b = b(:);
        for it = 1:40
            q = q(:); F = evalF(q, l); F = F(:);
            Nc  = tq(:).'*(q - qc) + tl*(l - lc) - ds;
            res = sqrt(F'*F + Nc^2);
            if res < arctol, conv = true; break; end
            [a, id] = mumps_bsolve(id, Jr, -F); a = a(:);
            dl = (-Nc - tq(:).'*a)/(tl + tq(:).'*b); dl = dl(1);
            dq = a + dl*b; q = q(:) + dq(:); l = l + dl;
        end
    end

    dt = toc(t_step)/60; cum = cum + dt;
    urms = sqrt(2)*l; amp = norm(q(1:N) + 1i*q(N+1:twoN));
    fprintf(' STEP %d: A0=%.6e u_rms=%.5e |q|=%.4e it=%d fac=%d age=%d %s (ds=%.3e) [%.1f min, cum %.1f, nfac_tot=%d]\n', ...
            step, l, urms, amp, it, nfac_step, age, tern(conv,'ok','FAIL'), ds, dt, cum, nfac_total);

    if ~conv
        if ~isempty(id), id = mumps_destroy(id); id = []; end
        ds = ds/2;
        if ds < 1e-6, fprintf(' ds too small; stop\n'); info.status = 'ds_underflow'; break; end
        continue
    end
    qc = q; lc = l; it_prev = it; age = age + 1;
    traj(end+1,:) = [l, urms, amp, it, nfac_step]; %#ok<AGROW>
    if ~isempty(ckfile)
        save(ckfile, 'q', 'l', 'N', 'twoN', 'step', 'traj', 'nfac_total', '-v7.3');
    end

    if l >= A_target - 1e-9
        if ~isempty(id), id = mumps_destroy(id); id = []; end
        fprintf('\n >>> reached target A0=%.6e u_rms=%.5f |q|=%.4e\n', l, sqrt(2)*l, amp);
        fprintf(' >>> TOTAL wall time: %.1f min (%.2f h) over %d steps | nfac_total=%d\n', ...
                toc(t_all)/60, toc(t_all)/3600, step, nfac_total);
        info.status = 'converged';
        break
    end

    tq_old = tq; tl_old = tl; ds_prev = ds;
    w = b; tln = 1/sqrt(w'*w + 1); tqn = tln*w;
    if (tqn'*tq + tln*tl) < 0, tqn = -tqn; tln = -tln; end
    tq = tqn; tl = tln;
    if it <= 6, ds = min(ds*1.2, dsmax); elseif it >= reuse_maxit-2, ds = ds*0.6; end
end
if ~isempty(id), id = mumps_destroy(id); end

info.nfac_total = nfac_total;
info.steps      = size(traj,1);
info.traj       = traj;
info.wall_min   = toc(t_all)/60;
info.lambda     = lc;

%% Unpack final state into StabRes
StabRes = SR0; Mode = Mode0;
q_complex = qc(1:N) + 1i*qc(N+1:twoN);
for k = 1:ctx0.nphys
    j   = ctx0.phys_modes(k);
    off = (k-1)*ctx0.bsize;
    qj  = q_complex(off+1:off+ctx0.bsize).';
    StabRes.q(j, :) = qj;
    jc = ctx0.nf + 1 - j;
    if jc ~= j; StabRes.q(jc, :) = conj(qj); end
end
[StabRes, ~, ~] = update_modes(SG0, Bnd0, Mode, Opt0, Stab0, StabRes);

% Apply the SAME postprocessing as the Picard path (main.m) so both solvers
% return an identical StabRes convention: shape functions normalised to unit
% max|u'|, amplitudes A = 2|u'|_max (mean-flow distortion = |u'|_max), and
% conjugate modes dropped. Uses the internal, fully-active Mode0 (RunJ=1:nf).
% main.m skips its own postprocess on the Newton path to avoid doubling this.
Init_pp.iu = 1i;
StabRes = postprocess(SG0, Init_pp, Mode, Opt0, StabRes);

fprintf('\n=== solve_newton DONE. %.1f min (%.2f h). nfac_total=%d, steps=%d, status=%s ===\n', ...
        info.wall_min, info.wall_min/60, nfac_total, info.steps, info.status);
end

%% Local helpers
function v = get_o(Opt, name, dflt)
if isfield(Opt, name) && ~isempty(Opt.(name)); v = Opt.(name); else; v = dflt; end
end
function s = tern(c, a, b)
if c, s = a; else, s = b; end
end

%% MUMPS helpers
%  JOB=4 (analyse+factor), JOB=3 (solve), JOB=-2 (destroy). The factor is held
%  across steps, so the id/Jr pair is threaded explicitly through the caller
%  rather than captured in a closure: refactoring destroys the instance, and a
%  closure over a destroyed id segfaults on the next solve.
function id = mumps_factor(Jr, Opt)
% Guard: a NaN/Inf in the Jacobian (e.g. the state overshot into a nonphysical
% rho'/T' excursion) makes dmumps hard-crash. Fail with a clear message instead.
% Check only stored nonzeros (structural zeros are finite; ~isfinite on the full
% sparse matrix would densify it and blow up memory).
nbad = nnz(~isfinite(nonzeros(Jr)));
if nbad > 0
    error('mumps_factor:badJacobian', ...
          ['Jacobian has %d non-finite entries before factorisation — the state ' ...
           'is nonphysical (likely a cold-start overshoot). Lower Opt.newton_start ' ...
           'or newton_dsmax.'], nbad);
end
% Attempt schedule, columns = [ICNTL(14) relax, ICNTL(22) out-of-core].
% In-core first, out-of-core last: the factor grows ~6x along the branch
% (21.9 GB at the ILST seed -> 123 GB at u_rms=4e-3, ny=40), so most steps fit
% in RAM and only the few near saturation pay the ~1.7x disk penalty. OOC drops
% the in-core need to roughly one front, so it cannot fail on memory; it spills
% to /tmp. Nothing else survived measurement: BLR does not compress this
% operator (13% mem, 1.5x slower, solve destroyed) and block/polynomial-
% preconditioned Krylov diverges at saturation. See project_m4_memory_wall.
sched = [100 0; 200 0; 100 1];
if isfield(Opt,'newton_mumps_relax') && ~isempty(Opt.newton_mumps_relax)
    r = Opt.newton_mumps_relax(:);
    sched = [[r, zeros(numel(r),1)]; [r(1), 1]];   % keep an OOC last resort
end
if isfield(Opt,'mumps_ooc') && ~isempty(Opt.mumps_ooc)
    sched(:,2) = Opt.mumps_ooc;                    % force OOC on/off throughout
end
info1 = NaN;
fprintf('   [mumps_factor] Jr: %dx%d nnz=%.3e square=%d sympat=%d maxabs=%.3e\n', ...
        size(Jr,1), size(Jr,2), nnz(Jr), size(Jr,1)==size(Jr,2), ...
        nnz(Jr-Jr.')==0, full(max(abs(nonzeros(Jr))))+0);
for attempt = 1:size(sched,1)
    relax = sched(attempt,1); ooc = sched(attempt,2);
    crashed = false; info2 = NaN;
    try
        id = initmumps; id.SYM = 0; id = dmumps(id);
        % ICNTL(7)=0 (AMD), not METIS: the saturated M=4 Jacobian does not fit
        % under METIS. Measured Jul 2026 on lsl, ny=40, A0=2.834e-3 (u_rms=4e-3),
        % real JOB=4:   AMD 9.29 min / 123 GB  vs  METIS 14.52 min / 216 GB.
        % The box has 188 GB, so METIS sits just over the line -> intermittent
        % INFOG(1)=-13. AMD leaves 65 GB headroom.
        % NB the seed (Stab.IC='ILST') Jacobian is ~6x sparser than the saturated
        % one (nnz 9.1e7 vs 5.5e8): ILST populates only the fundamental, so the
        % convolution blocks are empty until Newton fills them in on step 1.
        % Never size this solver from a seed-state probe. ICNTL(7)=4 (PORD)
        % segfaults this MUMPS 5.3.5 build.
        vrb = get_o(Opt,'newton_mumps_verbose',0);   % 0=silent; >0 -> MUMPS diagnostics
        id.ICNTL(3) = tern(vrb>0,6,0); id.ICNTL(4) = vrb; id.ICNTL(16) = 1;
        id.ICNTL(7) = get_o(Opt, 'newton_mumps_ordering', 0);
        id.ICNTL(14) = relax;
        id.ICNTL(22) = ooc;
        id.JOB = 4; id = dmumps(id, Jr);
    catch ME
        crashed = true;
        fprintf('   [mumps_factor] exception: %s | %s\n', ME.identifier, ME.message);
    end
    if ~crashed && isstruct(id) && ~isempty(id.INFOG)
        info1 = id.INFOG(1); info2 = id.INFOG(2);
        if info1 >= 0
            if attempt > 1
                fprintf('   [mumps_factor] OK on attempt %d (relax=%d ooc=%d)\n', attempt, relax, ooc);
            end
            return
        end
        fprintf('   [mumps_factor] attempt %d relax=%d ooc=%d FAILED: INFOG(1)=%d INFOG(2)=%d\n', ...
                attempt, relax, ooc, info1, info2);
        try, id.JOB = -2; id = dmumps(id); catch, end %#ok<CTCH>
    else
        fprintf('   [mumps_factor] attempt %d relax=%d ooc=%d HARD CRASH\n', attempt, relax, ooc);
        try, if isstruct(id), id.JOB = -2; id = dmumps(id); end; catch, end %#ok<CTCH>
    end
    pause(20);
end
error('mumps_factor: failed after %d attempts (last INFOG(1)=%g)', size(sched,1), info1);
end

function [X, id] = mumps_bsolve(id, Jr, RHS)
id.JOB = 3; id.RHS = RHS; id = dmumps(id, Jr); X = double(id.SOL);
end

function id = mumps_destroy(id)
id.JOB = -2; id = dmumps(id);
end
