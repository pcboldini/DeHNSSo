function [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS)
%MAIN  DeHNSSo main solver loop.
%
%   [StabRes, StabGrid]      = main(Stab, StabGrid, Opt)
%   [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt)       % also return LHS
%   [StabRes, StabGrid, LHS] = main(Stab, StabGrid, Opt, LHS)  % reuse LHS (Opt.rerun='on')
%
%   Solves the nonlinear Harmonic Navier-Stokes equations iteratively:
%   1. Initialise grid, operators, inflow, and boundary conditions
%   2. Convergence loop: solve modes, update amplitudes, check convergence
%   3. Post-process: normalize shape functions, compute wavenumbers
%
%   LHS struct fields: MLA1..MLA5, MFD_diff and LU_cache: base operators and
%   their per-mode factorisations. Pass LHS back in with Opt.rerun='on' to skip
%   both assembly and factorisation (omega/beta/A0 sweeps on a fixed base flow).

%% Input validation
if ~exist('StabGrid','var') || ~isstruct(StabGrid); StabGrid = struct(); end
if ~exist('Opt','var')      || ~isstruct(Opt);      Opt = struct();      end
if ~exist('Stab','var')     || ~isstruct(Stab);     Stab = struct();     end
if ~exist('LHS','var');     LHS = struct(); end
if ~isfield(Opt,'rerun');   Opt.rerun = 'off'; end
if ~isfield(Opt,'compressible'); Opt.compressible = 'off'; end
is_compr = strcmpi(Opt.compressible, 'on');

% Crash-safe cleanup of out-of-core block-Thomas factor files: sweep any streamed
% factors registered during THIS run on any in-MATLAB exit (error, Ctrl-C, normal
% return). Skipped when LHS is requested (nargout>=3), so a rerun keeps its factor.
bt_mark = block_thomas('free', [], 'mark');
if nargout < 3
    bt_ooc_cleanup = onCleanup(@() block_thomas('free', bt_mark, 'sweep')); %#ok<NASGU>
end

fprintf('\n');
fprintf(' ______     _   _  _   _  _____ _____       \n');
fprintf(' |  _  \\   | | | || \\ | |/  ___/  ___|      \n');
fprintf(' | | | |___| |_| ||  \\| |\\ `--.\\ `--.  ___  \n');
fprintf(' | | | / _ \\  _  || . ` | `--. \\`--. \\/ _ \\ \n');
fprintf(' | |/ /  __/ | | || |\\  |/\\__/ /\\__/ / (_) |\n');
fprintf(' |___/ \\___\\_| |_/\\_| \\_/\\____/\\____/ \\___/ \n');
fprintf('                                       v2    \n');
fprintf('  Delft Harmonic Navier-Stokes Solver\n');
fprintf('  =========================================\n');

%% Extract base-flow struct from StabGrid
fprintf('Extracting base flow ... ');
BF.Re  = StabGrid.Re;
BF.U   = StabGrid.U;    BF.V   = StabGrid.V;    BF.W   = StabGrid.W;
BF.dxU = StabGrid.dxU;  BF.dxV = StabGrid.dxV;  BF.dxW = StabGrid.dxW;
BF.dyU = StabGrid.dyU;  BF.dyV = StabGrid.dyV;  BF.dyW = StabGrid.dyW;

if is_compr
    required = {'rho','T','dxrho','dyrho','dxT','dyT', ...
                'dxxU','dxxV','dxxW','dxxT', ...
                'ddyU','ddyV','ddyW','ddyT','dxyU','dxyV', ...
                'Ma','Pr','gamma','Sstar','mu','dmuT','d2muT2'};
    missing = required(~cellfun(@(f) isfield(StabGrid,f), required));
    if ~isempty(missing)
        error('main:compressibleMissingFields', ...
            ['Opt.compressible=''on'' but StabGrid is missing: %s. ' ...
             'Re-run gridgen on a BF that carries rho and T.'], ...
            strjoin(missing, ', '));
    end
    for k = 1:numel(required)
        f = required{k};
        BF.(f) = StabGrid.(f);
    end
    fprintf('done (compressible, Re = %.2f, Ma = %.2f)\n', BF.Re, BF.Ma);
else
    fprintf('done (incompressible, 9 fields, Re = %.2f)\n', BF.Re);
end

%% 1. Initialisation
% solve_newton runs its own init/inflow (twice, to difference F w.r.t. the
% inflow amplitude), so it needs the RAW inputs, keep them before init overwrites.
StabGrid_raw = StabGrid; Stab_raw = Stab; Opt_raw = Opt;

[StabGrid, Init, Mode, Opt, Stab, StabRes] = init(StabGrid, Opt, Stab);
[Mode, StabRes] = inflow(Mode, Stab, StabGrid, StabRes, Opt);
[Bound] = setup_boundary(StabGrid, Opt, Stab);

%% Print configuration (after init sets defaults)
print_config(StabGrid, Stab, Opt);

%% 2. Convergence loop
% Reuse from a previous run (Opt.rerun='on'): restore the assembled blocks and
% the per-mode LU cache (both A0-independent) to skip assembly + factorisation.
MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MLA5 = []; MFD_diff = [];
if strcmpi(Opt.rerun,'on')
    has_mla = isfield(LHS,'MLA1') && ~isempty(LHS.MLA1);
    has_lu  = isfield(LHS,'LU_cache') && ~isempty(fieldnames(LHS.LU_cache));
    if has_mla
        MLA1     = LHS.MLA1;
        MLA2     = LHS.MLA2;
        if isfield(LHS, 'MLA3'); MLA3 = LHS.MLA3; end
        if isfield(LHS, 'MLA4'); MLA4 = LHS.MLA4; end
        if isfield(LHS, 'MLA5'); MLA5 = LHS.MLA5; end
        MFD_diff = LHS.MFD_diff;
        Opt.mmry = 1;   % signal solve_modes to skip LHS construction
    end
    if has_lu
        Opt.LU_cache = LHS.LU_cache;
        if isfield(LHS,'lu_order');       Opt.lu_order       = LHS.lu_order;       end
        if isfield(LHS,'lu_order_label'); Opt.lu_order_label = LHS.lu_order_label; end
    end
    if has_mla || has_lu
        if has_lu; nlu = numel(fieldnames(LHS.LU_cache)); else; nlu = 0; end
        fprintf(' Reusing LHS matrices from previous run (%d cached LU factorisation(s) reused).\n', nlu);
    end
end

is_linear = isfield(Opt,'linear') && strcmpi(Opt.linear,'on');
% Opt.solver is either 'picard' or 'newton' (= pseudo-arclength continuation,
% src/newton/solve_newton.m). There is no fixed-amplitude Newton: a cold-start
% solve does not converge at the amplitudes this solver targets, which is what the
% continuation exists for. Validate rather than fall through to Picard silently.
if isfield(Opt,'solver') && ~ismember(lower(Opt.solver), {'picard','newton'})
    error('main:badSolver', ...
        ['Opt.solver must be ''picard'' or ''newton'' (= pseudo-arclength ' ...
         'continuation, set Opt.newton_target); got ''%s''.'], Opt.solver);
end
is_newton = isfield(Opt,'solver') && strcmpi(Opt.solver, 'newton');

if is_compr && isfield(Opt,'adjoint') && strcmpi(Opt.adjoint,'on')
    error('main:compressibleAdjoint', ...
        'Adjoint not yet implemented for the compressible path.');
end

if is_newton
    %% Newton solve: pseudo-arclength continuation in A0_fund (see solve_newton.m)
    if ~isfield(Opt,'newton_target') || isempty(Opt.newton_target)
        error('main:noNewtonTarget', ...
              ['Opt.solver=''newton'' runs pseudo-arclength continuation and needs ' ...
               'Opt.newton_target (target A0_fund). See docs/user_guide.md.']);
    end
    fprintf('\n--- Newton solve (pseudo-arclength continuation) ---\n');
    t_loop = tic;
    [StabRes, info] = solve_newton(Stab_raw, StabGrid_raw, Opt_raw);
    fprintf('--- Newton %s after %d steps, %d factorisations (%.1f s) ---\n\n', ...
            info.status, info.steps, info.nfac_total, toc(t_loop));
    StabRes.newton_info = info;
    MLA1 = []; MLA2 = []; MLA3 = []; MLA4 = []; MLA5 = []; MFD_diff = [];
elseif is_linear
    %% Linear solve (single iteration, no NLT, no amplitude ramping)
    fprintf('\n--- Linear solve ---\n');
    t_loop = tic;
    Opt.dal = zeros(Mode.nf, 1);
    Opt.AF  = 1;

    [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff] = solve_modes( ...
        BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);
    [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes);

    if strcmpi(Opt.plot,'on')
        plot_intermediate(StabGrid, StabRes, Mode, Stab, Opt);
    end

    fprintf('--- Linear solve done (%.1f s) ---\n\n', toc(t_loop));
else
    %% Nonlinear convergence loop
    % AF-continuation Picard with damping and optional Anderson acceleration
    % (no Newton takeover: high amplitudes use Opt.solver='newton' instead)
    fprintf('\n--- Convergence loop ---\n');
    iter = 0;
    t_loop = tic;
    picard_stall_count = 0;
    AF_stall_check = -1;
    Opt.AF_history = [];
    iters_at_AF1 = 0;          % Picard iters spent at AF=1
    max_A_prev   = 0;          % previous iter's max_A — for runaway detection
    stall_msg_printed = false; % only print stall-disabled message once

    % Anderson acceleration state (Walker-Ni 2011). Set Opt.anderson = m to
    % enable (history depth, typical 3-6); 0 disables.
    if ~isfield(Opt,'anderson'); Opt.anderson = 0; end
    AA_on     = Opt.anderson > 0;
    AA_m      = Opt.anderson;
    AA_Xhist  = [];
    AA_Fhist  = [];
    AA_AF_ref = [];   % AF at which history was recorded — reset on AF change

    % Warm-start (Stab.IC='LOAD'): seed the NLT forcing from the loaded converged
    % field so the FIRST solve reproduces its saturated amplitude. Without this,
    % Stab.f = 0 on iter 1 -> a linear solve, and since under-relaxation only kicks
    % in at k1>0 the loaded field is overwritten, re-ramping the amplitude from A0.
    if strcmpi(Stab.IC,'LOAD') && strcmpi(Opt.allmodes,'on')
        [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes);
    end

    while max(abs(Opt.dal)) >= Opt.Conv || Opt.AF ~= 1
        iter = iter + 1;
        Opt.dal = zeros(Mode.nf, 1);
        t_iter = tic;

        % Damped Picard iteration
        q_old = StabRes.q;

        % Full Picard step: q_picard = ML⁻¹ * (rBC + f_NLT)
        [StabRes, Opt, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff] = solve_modes( ...
            BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff);

        % AVAL rescales the first iterate down to AMAX based on the
        % linear-regime prediction. Skip it when LOAD has warm-started
        % the full domain: the warm-start is already at physical saturation
        % and scaling it down throws away the continuation benefit.
        loaded_full_domain = strcmpi(Stab.IC,'LOAD');
        if nnz(Stab.A0) ~= 1 && Opt.k1 == 0 && ~loaded_full_domain
            [Opt, StabRes] = AVAL(StabGrid, Opt, Mode, StabRes);
        elseif loaded_full_domain && Opt.k1 == 0
            % With warm-start, default AF=1 from init.m stands; no damping.
            Opt.AVAL = max(abs(StabRes.q(:)));
        end

        % Anderson acceleration (optional): least-squares fixed-point
        % accelerator over the last m iterates. Falls back to plain Picard
        % if the accelerated step is worse.
        if AA_on && Opt.k1 > 0
            q_picard  = StabRes.q;                  % g(q_old)
            f_k_vec   = q_picard(:) - q_old(:);     % residual of fixed-point

            % Reset history if AF changed since last iterate
            if ~isempty(AA_AF_ref) && AA_AF_ref ~= Opt.AF
                AA_Xhist = []; AA_Fhist = [];
            end

            % Try Anderson step if we have at least 1 previous iterate
            if ~isempty(AA_Xhist)
                dX = [AA_Xhist(:, 2:end), q_old(:)]   - [AA_Xhist(:, 1:end)];
                dF = [AA_Fhist(:, 2:end), f_k_vec]    - [AA_Fhist(:, 1:end)];
                % Tikhonov-regularized LS: (dFᴴdF + λI)γ = dFᴴ f_k.
                % λ scales with ||dF||²_F so it's dimensionally consistent
                % and degrades gracefully to plain Picard when dF is
                % rank-deficient (parallel residuals near saturation).
                lambda = 1e-10 * norm(dF, 'fro')^2;
                gamma  = (dF'*dF + lambda * eye(size(dF,2))) \ (dF' * f_k_vec);
                q_AA_vec = q_picard(:) - (dX + dF) * gamma;

                % Safety: use AA only if result is finite AND not a
                % bigger blowup than plain Picard. Otherwise keep q_picard.
                max_AA = max(abs(q_AA_vec));
                max_P  = max(abs(q_picard(:)));
                if all(isfinite(q_AA_vec)) && max_AA < 1.5 * max_P
                    StabRes.q = reshape(q_AA_vec, size(q_picard));
                    fprintf('   [Anderson: m=%d, γ=(%s)]\n', ...
                        size(dX,2), sprintf('%.2f ', gamma));
                end
            end

            % Push new entry, drop oldest if full
            AA_Xhist(:, end+1) = q_old(:);        %#ok<AGROW>
            AA_Fhist(:, end+1) = f_k_vec;         %#ok<AGROW>
            if size(AA_Xhist, 2) > AA_m
                AA_Xhist(:, 1) = [];
                AA_Fhist(:, 1) = [];
            end
            AA_AF_ref = Opt.AF;
        end

        % Damp if amplitudes grow too fast (prevents blow-up at high AF).
        % Applied AFTER Anderson, as a final safety net.
        if Opt.k1 > 0
            q_picard = StabRes.q;
            dq = q_picard - q_old;
            max_A_old = max(max(abs(q_old(Mode.RunJ, :))));
            max_A_new = max(max(abs(q_picard(Mode.RunJ, :))));

            if max_A_new > 2 * max(max_A_old, 1e-6)
                alpha = max_A_old / max(max_A_new, eps);
                alpha = max(alpha, 0.1);
                StabRes.q = q_old + alpha * dq;
                fprintf('   [damped Picard: alpha=%.2f]\n', alpha);
            end
        end

        % Fixed under-relaxation (Opt.relax = omega in (0,1)): q <- q_old + omega*(q - q_old).
        % Reduces the effective spectral radius so a marginally-repelling Picard fixed
        % point becomes attracting. Use with Anderson off (they fight each other).
        if Opt.k1 > 0 && isfield(Opt,'relax') && ~isempty(Opt.relax) && Opt.relax < 1
            StabRes.q = q_old + Opt.relax * (StabRes.q - q_old);
        end

        % Update modes, amplitudes, NLT forcing
        [StabRes, Mode, Stab] = update_modes(StabGrid, Bound, Mode, Opt, Stab, StabRes);

        % Check convergence and ramp AF
        AF_before = Opt.AF;
        [Opt, Stab, StabRes] = check_convergence(StabGrid, Bound, Mode, Opt, Stab, StabRes);
        Opt.AF_history(end+1) = Opt.AF;

        % Status
        n_active = numel(Mode.RunJ(Mode.RunJ >= round(Mode.nf/2)));
        max_A = max(StabRes.A(:));
        max_dal = max(abs(Opt.dal));
        solver_tag = 'P';
        fprintf('  Iter %3d [%s] | AF=%.3f | %d modes | max(A)=%.3e | max(dal)=%.3e | %.1f s\n', ...
            iter, solver_tag, Opt.AF, n_active, max_A, max_dal, toc(t_iter));

        % --- Picard stall detection ---
        % Reset counter when AF changes or number of active modes changes
        if ~isfield(Opt, 'n_active_prev'); Opt.n_active_prev = n_active; end
        modes_changed = (n_active ~= Opt.n_active_prev);
        Opt.n_active_prev = n_active;

        if modes_changed
            picard_stall_count = 0;
            picard_af_count = 0;
            AF_stall_check = Opt.AF;
            dal_prev = inf;
        elseif Opt.AF ~= AF_before
            % AF just ramped, reset
            picard_stall_count = 0;
            picard_af_count = 0;
            AF_stall_check = Opt.AF;
            dal_prev = inf;
        elseif max_dal > Opt.Conv
            if ~exist('dal_prev','var'); dal_prev = inf; end
            if ~exist('picard_af_count','var'); picard_af_count = 0; end
            picard_af_count = picard_af_count + 1;  % total iters at this AF
            if max_dal >= dal_prev
                picard_stall_count = picard_stall_count + 1;
            else
                picard_stall_count = 0;
            end
            dal_prev = max_dal;
        else
            picard_stall_count = 0;
            picard_af_count = 0;
            dal_prev = max_dal;
        end

        % Track time at AF=1 and detect runaway early, regardless of dal
        % (which gets masked to 0 when AVAL >= 2). These catch Picard
        % *before* it loses the physical saturation state.
        if Opt.AF >= 1 - eps
            iters_at_AF1 = iters_at_AF1 + 1;
        else
            iters_at_AF1 = 0;
        end
        runaway_jump = (iters_at_AF1 > 1) && (max_A > 1.5 * max_A_prev) && (max_A > 0.15);
        max_A_prev = max_A;

        % Divergence, NOT a stall, and not recoverable. max_A > 1 means the
        % perturbation exceeds the freestream (unphysical for this expansion);
        % NaN/Inf means it has already overflowed. Continuing from here lets the
        % amplitudes reset and dal collapse to exactly 0, which check_convergence
        % then reports as convergence, returning a plausible-looking field at
        % ~A0 that is garbage. Stop instead of faking success.
        if ~isfinite(max_A) || max_A > 1
            error('main:picardDiverged', ...
                ['Picard diverged at iteration %d (AF=%.3f): max(A) = %.3e. ' ...
                 'This is divergence, not a stall — the iteration cannot recover. ' ...
                 'Try Opt.anderson = 5 (often enough near saturation), a lower ' ...
                 'Stab.A0_fund, or Opt.solver = ''newton'' (pseudo-arclength ' ...
                 'continuation) to reach this amplitude.'], iter, Opt.AF, max_A);
        end

        % Stall (recoverable: the solution is still finite, Picard is just slow):
        %  - dal not decreasing for 5 iters, OR
        %  - too many iters at same AF (Picard too slow), OR
        %  - stall_count + max_A cross threshold, OR
        %  - Picard has been at AF=1 for 5+ iters AND dal has stopped
        %    decreasing for 3+ iters (done its AF-ramp job, truly stuck), OR
        %  - runaway signature: max_A jumps >50% in one iter at AF=1.
        if picard_stall_count >= 5 || picard_af_count >= 15 || ...
           (picard_stall_count >= 3 && max_A > 0.1) || ...
           (iters_at_AF1 >= 5 && picard_stall_count >= 3) || ...
           runaway_jump
            picard_stall_count = 10;
        end

        if picard_stall_count >= 10
            % No Newton takeover: the solver is either Picard or Newton
            % (Opt.solver='newton' = pseudo-arclength continuation), never a
            % mix. Report the stall once and keep running Picard.
            picard_stall_count = 0;
            if ~stall_msg_printed
                fprintf(['\n  ** Picard stalled at AF=%.3f — continuing Picard. ' ...
                         'For high amplitudes use Opt.solver=''newton'' ' ...
                         '(pseudo-arclength continuation) instead. **\n\n'], Opt.AF);
                stall_msg_printed = true;
            end
        end

        % Plotting
        if strcmpi(Opt.plot,'on')
            plot_intermediate(StabGrid, StabRes, Mode, Stab, Opt);
        end
        if Opt.Sweep == 1 && max(abs(Opt.dal)) <= Opt.Sweep * Opt.Conv
            [Opt, StabRes] = save_sweep(StabGrid, Mode, Opt, Stab, StabRes);
        end
        if strcmpi(Opt.plot,'on')
            plot_convergence(StabGrid, Mode, Stab, Opt);
        end

        % Per-iteration checkpoint: persist the latest finite iterate so a kill,
        % stall (no max-iter cap on this loop), or divergence-error mid-run does
        % not lose the whole run (the step-3 save only fires on a clean loop exit).
        % Only overwrite when q is finite, so a good checkpoint survives a blow-up.
        if isfield(Opt,'save_checkpoint') && ~isempty(Opt.save_checkpoint) ...
           && all(isfinite(StabRes.q(:)))
            ck_dir = fileparts(Opt.save_checkpoint);
            if ~isempty(ck_dir) && ~exist(ck_dir,'dir'); mkdir(ck_dir); end
            save(Opt.save_checkpoint, 'StabRes', 'StabGrid', 'Stab', 'Mode');
        end
    end

    % Report. Distinguish true convergence from the Newton-stall bail-out.
    converged = (max(abs(Opt.dal)) < Opt.Conv) && (Opt.AF == 1);
    if converged; label = 'Converged'; else; label = 'Stopped (did NOT converge)'; end
    % Expose Picard outcome to callers (mirrors StabRes.newton_info for the
    % Newton path). Bench/diagnostic scripts can check .picard_info.status.
    if converged
        picard_status = 'converged';
    else
        picard_status = 'stalled';
    end
    StabRes.picard_info = struct( ...
        'status',   picard_status, ...
        'iters',    iter, ...
        'max_dal',  max(abs(Opt.dal)), ...
        'AF_final', Opt.AF);
    try
        [~, rss] = system(sprintf('ps -o rss= -p %d', feature('getpid')));
        fprintf('--- %s in %d iterations (%.1f s, peak %.1f GB) ---\n\n', ...
            label, iter, toc(t_loop), str2double(strtrim(rss))/1e6);
    catch
        fprintf('--- %s in %d iterations (%.1f s) ---\n\n', label, iter, toc(t_loop));
    end
end

%% 3. Save checkpoint (optional, before postprocess)
% Dumps StabRes with full mode set + .q intact, suitable for Stab.IC='LOAD'
% in a subsequent run (amplitude continuation, parameter sweep). Triggered by
% Opt.save_checkpoint = 'path/to/file.mat'. Only saves if solution is finite.
if isfield(Opt,'save_checkpoint') && ~isempty(Opt.save_checkpoint) ...
   && all(isfinite(StabRes.q(:)))
    ckpt_path = Opt.save_checkpoint;
    ckpt_dir  = fileparts(ckpt_path);
    if ~isempty(ckpt_dir) && ~exist(ckpt_dir,'dir'); mkdir(ckpt_dir); end
    save(ckpt_path, 'StabRes', 'StabGrid', 'Stab', 'Mode');
    fprintf('--- Checkpoint saved: %s ---\n', ckpt_path);
end

%% 4. Post-processing
half = round(Mode.nf/2);
% Skip postprocess if the solve diverged (StabRes contains NaN/Inf). Running
% interp1 / spline fits on NaN would crash. Just return the raw (broken)
% StabRes so the user can inspect.
if is_newton
    % solve_newton already applied postprocess internally (with its own active
    % Mode); running it again here would double-normalise and re-drop conjugates.
elseif any(~isfinite(StabRes.q(:)))
    warning('main:postprocess_skipped', ...
        'Solve diverged (non-finite values in StabRes.q) — skipping postprocess.');
else
    StabRes = postprocess(StabGrid, Init, Mode, Opt, StabRes);
end

%% Save output (optional): final post-processed results for later analysis/plotting.
% Opt.save_output = 'path/to/file.mat' (distinct from Opt.save_checkpoint, which
% dumps the pre-postprocess state for Stab.IC='LOAD' continuation).
if isfield(Opt,'save_output') && ~isempty(Opt.save_output)
    out_path = Opt.save_output;
    out_dir  = fileparts(out_path);
    if ~isempty(out_dir) && ~exist(out_dir,'dir'); mkdir(out_dir); end
    save(out_path, 'StabRes', 'StabGrid', 'Stab');
end

%% Pack LHS output (for reuse in subsequent runs)
LHS.MLA1     = MLA1;
LHS.MLA2     = MLA2;
LHS.MLA3     = MLA3;
LHS.MLA4     = MLA4;
LHS.MLA5     = MLA5;
LHS.MFD_diff = MFD_diff;
% Carry the LU cache so the next Opt.rerun can skip factorisation (reference, no copy).
if isfield(Opt,'LU_cache') && ~isempty(fieldnames(Opt.LU_cache))
    LHS.LU_cache = Opt.LU_cache;
    if isfield(Opt,'lu_order');       LHS.lu_order       = Opt.lu_order;       end
    if isfield(Opt,'lu_order_label'); LHS.lu_order_label = Opt.lu_order_label; end
end
% Out-of-core factor files are freed by the onCleanup guard installed at entry
% (crash-safe; skipped for nargout>=3 so a rerun keeps its factor).

%% Completion message
if isfield(Opt,'save_output') && ~isempty(Opt.save_output)
    fprintf('\n=== DeHNSSo finished — output saved to %s ===\n', Opt.save_output);
else
    fprintf('\n=== DeHNSSo finished (output not saved; set Opt.save_output to save) ===\n');
end

end
