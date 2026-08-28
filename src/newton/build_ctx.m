function ctx = build_ctx(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes)
%BUILD_CTX  Assemble the Newton solver context.
%
%   ctx = build_ctx(BF, Bound, StabGrid, Init, Mode, Opt, Stab, StabRes)
%
%   Builds the per-mode LHS operators, their block-preconditioner LU factors and
%   the inhomogeneous boundary-condition vectors, packed into one context struct.
%   Used by solve_newton (pseudo-arclength continuation in A0_fund) and by the
%   build_JNL instrumentation in paper/checks/.

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
% Keep only the difference: the MFD (0,0) operator differs from MLA1 in its
% boundary rows alone, so it is rebuilt below as MLA1 + MFD_diff.
MFD_diff = MLAMFD - MLA1;
clear MLAMFD

fprintf(' %.1fs\n', toc(t));

%% Active modes and context fields
% Conjugate modes are stored in the lower half of the spectrum, so the physical
% ones are those at or above the midpoint.
half  = round(Mode.nf/2);
phys  = Mode.RunJ(Mode.RunJ >= half);   % active physical modes only
nphys = numel(phys);

% Sizes. bsize is the length of one mode's state, so the full Newton state is
% nphys*bsize complex, twice that in the real form the solver works in.
ctx.phys_modes = phys;
ctx.nphys      = nphys;
ctx.nf         = Mode.nf;
ctx.ny         = StabGrid.ny;
ctx.nx         = StabGrid.nx;
ctx.nvar       = nvar;
ctx.bsize      = nvar*ctx.ny*ctx.nx;

% Grid and wall-normal operator, read by the NLT and by build_JNL.
ctx.StabGrid   = StabGrid;
ctx.D1         = StabGrid.D1;

% HB is the mode-interaction table: HB(:,:,j) lists the (gh,jk) pairs whose
% product forces mode j, which is what build_JNL walks.
ctx.HB         = Stab.HB;
ctx.omegavec   = StabRes.omegavec;
ctx.betavec    = StabRes.betavec;

% Inputs carried through untouched, so the residual can be re-evaluated at any
% lambda without the caller passing them again.
ctx.Bound      = Bound;
ctx.Mode       = Mode;
ctx.BF_orig    = BF;
ctx.Init       = Init;
ctx.Stab       = Stab;
ctx.StabRes    = StabRes;
ctx.Opt        = Opt;

% Spectral indices, and the curvilinear metrics used to form physical-space
% derivatives from the (xi,eta) operators.
ctx.Mvec       = Stab.Mvec;
ctx.Nvec       = Stab.Nvec;
ctx.dxi        = StabGrid.x(2) - StabGrid.x(1);   % legacy scalar, for display only
ctx.D_xi       = StabGrid.D_xi;                    % sparse non-uniform FD matrix in xi
ctx.xix        = StabGrid.xix;
ctx.etax       = StabGrid.etax;
ctx.xiy        = StabGrid.xiy;
ctx.etay       = StabGrid.etay;

% Per-mode operator and boundary-condition right-hand side, filled below.
ctx.ML  = cell(nphys, 1);
ctx.rBC = cell(nphys, 1);

Opt_rhs = Opt;
Opt_rhs.AF = 1;   % store UNIT amplitude BC; residual_hns scales by ctx.AF
key_cache = struct();

%% Per-mode LHS blocks
% The LHS is polynomial in (omega, beta), so each mode's operator is combined from
% the precomputed blocks rather than reassembled. Modes sharing an (omega, beta)
% pair share one operator through key_cache.
for k = 1:nphys
    j     = phys(k);
    omega = StabRes.omegavec(j);
    beta  = StabRes.betavec(j);
    key   = ml_key(omega, beta);

    if ~isfield(key_cache, key)
        % MLA5 is the i*beta*d_x cross-viscous block; compressible only, empty otherwise.
        b5 = 0; if is_compr && ~isempty(MLA5); b5 = beta; end
        % The branches are the one polynomial with its zero terms dropped, so a 2D
        % or steady mode does not pay for sparse additions that vanish anyway.
        if omega == 0 && beta == 0
            ML = MLA1 + MFD_diff;              % (0,0) MFD: its own boundary rows
        elseif omega == 0
            ML = MLA1 + MLA3*beta + MLA4*beta^2;
            if b5, ML = ML + MLA5*b5; end
        elseif beta == 0
            ML = MLA1 + MLA2*omega;
        else
            ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2;
            if b5, ML = ML + MLA5*b5; end
        end
        key_cache.(key) = struct('ML', ML);
    end
    ctx.ML{k} = key_cache.(key).ML;
    % Inflow BC at unit amplitude; residual_hns scales it by the continuation AF.
    ctx.rBC{k} = RHS(Bound, StabGrid, Mode, Opt_rhs, StabRes, j);
end

end


function key = ml_key(omega, beta)
% Cache key for one (omega, beta) pair; the (0,0) mean-flow distortion is named
% separately since it carries its own boundary rows.
if omega == 0 && beta == 0
    key = 'MFD';
else
    key = matlab.lang.makeValidName(sprintf('o%g_b%g', omega, beta));
end
end

function s = ternary_c(c,a,b); if c, s=a; else, s=b; end; end
