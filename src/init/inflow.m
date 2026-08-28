function [Mode,StabRes] = inflow(Mode,Stab,StabGrid,StabRes,Opt)
%INFLOW  Seed the inflow disturbance for all active modes.
%
%   [Mode, StabRes] = inflow(Mode, Stab, StabGrid, StabRes, Opt)
%
%   Builds the disturbance at the inlet station, one mode at a time. The
%   method is chosen by Stab.IC: 'ILST'/'CLST' (incompressible/compressible
%   LST, linear stability theory) compute a local eigenmode (via local
%   functions ilst_inflow/clst_inflow), 'PSE' (parabolized stability
%   equations) marches a non-parallel seed instead (via
%   solver_IPSE/solver_CPSE), 'LOAD' warm-starts from a saved file, and
%   'ZERO' leaves the inlet undisturbed so the field is generated entirely
%   by wall/freestream forcing. Each mode is scaled by Stab.A0, given a
%   phase reference (wall pressure or max|u|), and written into
%   StabRes.u/v/w/p(/rho/T)/phi at station 1. Mode.RunJ is then updated to
%   list the active modes (those seeded here plus those forced at the
%   wall).

nonzeromodeswall = [];
nonzeromodesic   = [];

% Modes activated by wall/freestream actuation (Stab.bcw / Stab.bct)
nf_bc = min(size(Stab.bcw,3), Mode.nf);
for j = 1:nf_bc
    if any(Stab.bcw(:,:,j),'all') || any(Stab.bct(:,:,j),'all')
        nonzeromodeswall(end+1) = j;
    end
end

is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
nvar     = 4 + double(is_compr);

% Amplitude convention: Stab.A0 (and Stab.bcw, halved in RHS) is the PHYSICAL peak
% amplitude 2|û| of the real disturbance q' = q̂ e^{-iωt} + c.c. The inflow seeds the
% complex amplitude û = A0/2; post-process then reports StabRes.A = 2|û| = A0, so the
% amplitude you set in equals the amplitude you get out (see docs/user_guide.md §5).
A0s = 0.5 * Stab.A0;

StabRes.phi    = zeros(Mode.nf,nvar*StabGrid.ny,StabGrid.nx);
StabRes.alpha  = zeros(Mode.nf,StabGrid.nx);
StabRes.A      = zeros(Mode.nf,StabGrid.nx);
StabRes.u      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.v      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.w      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.p      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
if is_compr
    StabRes.rho = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
    StabRes.T   = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
end

if is_compr && ~(strcmpi(Stab.IC,'CLST') || strcmpi(Stab.IC,'LOAD') || strcmpi(Stab.IC,'PSE') || strcmpi(Stab.IC,'ZERO'))
    error('inflow:badIC', ...
        'Opt.compressible=''on'' requires Stab.IC=''CLST'', ''PSE'', ''LOAD'' or ''ZERO''; got ''%s''.', Stab.IC);
end
if ~is_compr && strcmpi(Stab.IC,'CLST')
    error('inflow:badIC', 'Stab.IC=''CLST'' requires Opt.compressible=''on''.');
end

% Freestream BC flag for the LST inflow solvers: single 'H_DR'/'H_NM', Neumann
% if any component is. clst_inflow uses it; the ILST seed decays so it ignores it.
if ~exist('Opt','var') || ~isfield(Opt,'bc_top')
    bc_top = 'H_DR';
elseif iscell(Opt.bc_top)
    if any(strcmpi(Opt.bc_top, 'H_NM'))
        bc_top = 'H_NM';
    else
        bc_top = 'H_DR';
    end
else
    bc_top = Opt.bc_top;
end

if Stab.IC=="LOAD"
    % phiIC width = nvar*ny: 4 incompressible (u,v,w,p), 5 compressible (u,v,w,ρ,T).

    assert(isfield(Stab,'ICfile') && ~isempty(Stab.ICfile), ...
        'Stab.ICfile must be set when Stab.IC = ''LOAD''.');
    assert(isfile(Stab.ICfile), 'ICfile not found: %s', Stab.ICfile);

    loaded = load(Stab.ICfile);

    if isfield(loaded,'phiIC')
        phiIC_load = loaded.phiIC;
        ny = StabGrid.ny;

        % Per-variable layout, the natural one for an externally computed
        % profile: [ny x nvar] for a single mode, or [ny x nvar x k] for k of
        % them, with the wall-normal index down each column and one column per
        % variable (u,v,w,p or u,v,w,rho,T). Folded here into the internal
        % [k x nvar*ny] row stacking; MATLAB is column-major, so A(:) already
        % concatenates the columns in the required order.
        if size(phiIC_load,1) == ny && size(phiIC_load,2) == nvar
            nblk = size(phiIC_load,3);
            stacked = zeros(nblk, nvar*ny);
            for r = 1:nblk
                A = phiIC_load(:,:,r);
                stacked(r,:) = A(:).';
            end
            phiIC_load = stacked;
        end

        assert(size(phiIC_load,2) == nvar*StabGrid.ny, ...
            'Loaded phiIC has %d columns but %d are expected (nvar*ny = %d*%d).', ...
            size(phiIC_load,2), nvar*StabGrid.ny, nvar, StabGrid.ny);

        % Two accepted forms. Without phiIC_modes the array spans the whole mode
        % list, row j being mode j. With it, only the listed modes are supplied:
        % phiIC_modes is a k-by-2 list of (m,n) pairs naming each row, so a linear
        % run needs a single row rather than a mostly empty nf-row array.
        if isfield(loaded,'phiIC_modes')
            pm = loaded.phiIC_modes;
            assert(size(pm,2) == 2 && size(pm,1) == size(phiIC_load,1), ...
                ['phiIC_modes must be [%d x 2], one (m,n) pair per phiIC row; ' ...
                 'got [%d x %d].'], size(phiIC_load,1), size(pm,1), size(pm,2));
            assert(isfield(Stab,'Mvec') && isfield(Stab,'Nvec'), ...
                'phiIC_modes needs Stab.Mvec/Nvec (set by init).');
            rows = zeros(size(pm,1),1);
            for r = 1:size(pm,1)
                jr = find(Stab.Mvec(:).' == pm(r,1) & Stab.Nvec(:).' == pm(r,2));
                assert(isscalar(jr), ...
                    'phiIC_modes row %d names mode (%d,%d), which is not in this run (M=%d, N=%d).', ...
                    r, pm(r,1), pm(r,2), Stab.M, Stab.N);
                rows(r) = jr;
            end
            assert(numel(unique(rows)) == numel(rows), 'phiIC_modes names the same mode twice.');
            full = zeros(Mode.nf, nvar*ny);
            full(rows,:) = phiIC_load;
            phiIC_load = full;
            fprintf('  LOAD: %d user profile(s) mapped to mode(s)%s\n', numel(rows), ...
                    sprintf(' (%d,%d)', pm.'));
        else
            assert(size(phiIC_load,1) == Mode.nf, ...
                ['Loaded phiIC has %d rows but %d are expected (the full mode list). ' ...
                 'Supply phiIC_modes, a k-by-2 list of (m,n) pairs, to give only some modes.'], ...
                size(phiIC_load,1), Mode.nf);
        end

        for j = 1:Mode.nf
            StabRes.u(j,:,1) = phiIC_load(j,      1:ny);
            StabRes.v(j,:,1) = phiIC_load(j,   ny+1:2*ny);
            StabRes.w(j,:,1) = phiIC_load(j, 2*ny+1:3*ny);
            if is_compr
                StabRes.rho(j,:,1) = phiIC_load(j, 3*ny+1:4*ny);
                StabRes.T(j,:,1)   = phiIC_load(j, 4*ny+1:5*ny);
            else
                StabRes.p(j,:,1)   = phiIC_load(j, 3*ny+1:4*ny);
            end
        end

    elseif isfield(loaded,'StabRes')
        SR = loaded.StabRes;
        % ICfile must be a CHECKPOINT (Opt.save_checkpoint): .q intact, full mode
        % set, shapes unnormalised. A results file (Opt.save_output) has been
        % through postprocess (strips q/phi/R), so it cannot restore solver state.
        if ~isfield(SR,'q')
            error('inflow:ICfileNotCheckpoint', ...
                ['ICfile "%s" has no StabRes.q — it looks like an Opt.save_output ' ...
                 'results file (postprocess strips q/phi/R). Stab.IC=''LOAD'' needs a ' ...
                 'checkpoint: set Opt.save_checkpoint in the producing run and load that.'], ...
                 Stab.ICfile);
        end
        assert(size(SR.u,1) == Mode.nf && size(SR.u,2) == StabGrid.ny, ...
            ['Loaded StabRes.u size [%d x %d x ...] does not match expected [%d x %d x ...]. ' ...
             'The checkpoint must come from a run with the same M/N and grid.'], ...
            size(SR.u,1), size(SR.u,2), Mode.nf, StabGrid.ny);
        if is_compr
            assert(isfield(SR,'rho') && isfield(SR,'T'), ...
                'Loaded StabRes lacks .rho/.T fields (compressible LOAD requires them).');
        end
        % Continuation: rescale the warm-start linearly when Stab.A0_fund changed.
        scale = 1;
        if isfield(loaded,'Stab') && isfield(loaded.Stab,'A0_fund') ...
                && isfield(Stab,'A0_fund') && loaded.Stab.A0_fund > 0
            scale = Stab.A0_fund / loaded.Stab.A0_fund;
        end
        if is_compr; flds = {'u','v','w','rho','T'}; else; flds = {'u','v','w','p'}; end
        if size(SR.u,3) == StabGrid.nx
            for j = 1:Mode.nf
                for f = 1:numel(flds)
                    StabRes.(flds{f})(j,:,:) = scale * SR.(flds{f})(j,:,:);
                end
            end
            if isfield(SR,'q') && size(SR.q,2) == nvar*StabGrid.ny*StabGrid.nx
                StabRes.q = scale * SR.q;
            end
            fprintf('  LOAD: full-domain warm-start (%d modes -> %d incl. conjugates, %d x %d), A0 scale = %.3f\n', ...
                    size(SR.u,1), Mode.nf, size(SR.u,2), size(SR.u,3), scale);
        else
            for j = 1:Mode.nf
                for f = 1:numel(flds)
                    StabRes.(flds{f})(j,:,1) = SR.(flds{f})(j,:,1);
                end
            end
            fprintf('  LOAD: inflow-only (station 1) — no downstream warm-start\n');
        end

    else
        error('ICfile must contain a variable named ''StabRes'' or ''phiIC''.');
    end

    if is_compr
        % Derive p' from EOS:  p' = (ρ'·T̄ + ρ̄·T') / (γ·Ma²)
        gM2 = StabGrid.gamma * StabGrid.Ma^2;
        for j = 1:Mode.nf
            StabRes.p(j,:,1) = (StabRes.rho(j,:,1) .* StabGrid.T(:,1).' ...
                              + StabGrid.rho(:,1).' .* StabRes.T(j,:,1)) / gM2;
        end
    end

    [~, nonzeromodesic] = find(Stab.A0 > 0);
    for j = flip(nonzeromodesic)
        if j >= round(Mode.nf/2)
            % Normalise by max(|u|), then scale by A0.
            y_int  = linspace(StabGrid.y(1), StabGrid.y(end), 4000);
            u1int  = interp1(StabGrid.y, abs(StabRes.u(j,:,1)), y_int, 'spline');
            ampltd = max(abs(u1int));
            if ampltd > 0
                StabRes.u(j,:,1) = StabRes.u(j,:,1) / ampltd;
                StabRes.v(j,:,1) = StabRes.v(j,:,1) / ampltd;
                StabRes.w(j,:,1) = StabRes.w(j,:,1) / ampltd;
                if is_compr
                    StabRes.rho(j,:,1) = StabRes.rho(j,:,1) / ampltd;
                    StabRes.T(j,:,1)   = StabRes.T(j,:,1)   / ampltd;
                else
                    StabRes.p(j,:,1) = StabRes.p(j,:,1) / ampltd;
                end
            end
            StabRes.u(j,:,1) = A0s(j) .* StabRes.u(j,:,1);
            StabRes.v(j,:,1) = A0s(j) .* StabRes.v(j,:,1);
            StabRes.w(j,:,1) = A0s(j) .* StabRes.w(j,:,1);
            if is_compr
                StabRes.rho(j,:,1) = A0s(j) .* StabRes.rho(j,:,1);
                StabRes.T(j,:,1)   = A0s(j) .* StabRes.T(j,:,1);
                StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.rho(j,:,1) StabRes.T(j,:,1)];
            else
                StabRes.p(j,:,1) = A0s(j) .* StabRes.p(j,:,1);
                StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.p(j,:,1)];
            end
            switch lower(Stab.phaseRef)
                case 'pwall'
                    theta0 = angle(StabRes.p(j,end,1));
                case 'umax'
                    [~,iy] = max(abs(StabRes.u(j,:,1)));
                    theta0 = angle(StabRes.u(j,iy,1));
                otherwise
                    error('Unknown phaseRef option');
            end
            StabRes.phi(j,:,1) = StabRes.phi(j,:,1) * exp(-1i*theta0);
        else
            % Conjugate mode
            j2 = round(Mode.nf/2) - (j - round(Mode.nf/2));
            StabRes.u(j,:,1)   = conj(StabRes.u(j2,:,1));
            StabRes.v(j,:,1)   = conj(StabRes.v(j2,:,1));
            StabRes.w(j,:,1)   = conj(StabRes.w(j2,:,1));
            if is_compr
                StabRes.rho(j,:,1) = conj(StabRes.rho(j2,:,1));
                StabRes.T(j,:,1)   = conj(StabRes.T(j2,:,1));
            else
                StabRes.p(j,:,1)   = conj(StabRes.p(j2,:,1));
            end
            StabRes.alpha(j,1) = -conj(StabRes.alpha(j2,1));
            StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
        end
    end
    fprintf('Initialising with LOAD solution from: %s\n', Stab.ICfile);

elseif Stab.IC=="ILST"
    [~,nonzeromodesic]= find(Stab.A0>0);

    for j = flip(nonzeromodesic) % Compute the inflow condition for every nonzero mode
    if j >= round(Mode.nf/2)% Physical mode (not its complex conjugate)
        [StabRes]=ilst_inflow(j,1,StabGrid,StabRes,bc_top,Opt);

        % ilst_inflow returns the shape normalised to unit-max|u'|; scale it by the physical amplitude
        StabRes.u(j,:,1) = A0s(j).*StabRes.u(j,:,1);
        StabRes.v(j,:,1) = A0s(j).*StabRes.v(j,:,1);
        StabRes.w(j,:,1) = A0s(j).*StabRes.w(j,:,1);
        StabRes.p(j,:,1) = A0s(j).*StabRes.p(j,:,1);
        StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.p(j,:,1)];
        switch lower(Stab.phaseRef)  % phase reference: 'pwall' or 'umax'
        case 'pwall'
            theta0 = angle(StabRes.p(j,end,1));  % phase of the wall pressure
        case 'umax'
            [~,iy] = max(abs(StabRes.u(j,:,1))); % location of max |u'|
            theta0 = angle(StabRes.u(j,iy,1));
        otherwise
            error('Unknown phaseRef option');
        end
        % Rotate the inflow perturbation to this phase reference
        StabRes.phi(j,:,1) = StabRes.phi(j,:,1) * exp(-1i*theta0);

    else % Conjugate mode: build it by symmetry from its physical partner
        j2 = round(Mode.nf/2)-(j-round(Mode.nf/2));% Index of the physical partner mode

        StabRes.u(j,:,1) = conj(StabRes.u(j2,:,1));
        StabRes.v(j,:,1) = conj(StabRes.v(j2,:,1));
        StabRes.w(j,:,1) = conj(StabRes.w(j2,:,1));
        StabRes.p(j,:,1) = conj(StabRes.p(j2,:,1));
        StabRes.alpha(j,1) = -conj(StabRes.alpha(j,1)); % conjugate mode: alpha_r flips sign, alpha_i (growth rate) unchanged
        StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
    end
    end
    
    if max(abs(Stab.bcw(2:end,:)))==0
        fprintf('Initialising with ILST solution at inflow. \n')
    else
        fprintf('Initialising with ILST solution at inflow and non-homogeneous b.c. at wall. \n')
    end

elseif Stab.IC=="CLST" || Stab.IC=="PSE"
    [~,nonzeromodesic] = find(Stab.A0>0);
    sym_inflow = isfield(Opt,'sym_inflow') && strcmpi(Opt.sym_inflow,'on');
    is_pse = strcmpi(char(Stab.IC),'PSE');   % linear-PSE (non-parallel) seed vs local CLST
    ov_ = real(StabRes.omegavec(:)); bv_ = real(StabRes.betavec(:)); tol_ = 1e-9;

    if sym_inflow
        % Symmetry-correct oblique seeding: CLST-solve only the true fundamental
        % (omega>0 & beta>0); derive the other seeded modes as its images on the
        % 2D base flow (Wbar=0): spanwise mirror flips w', reality conjugation
        % conjugates with alpha->-conj(alpha). Lands the (-1,1) partner that the
        % 1D index-mirror left at zero.
        is_canon = (ov_ >= -tol_) & (bv_ >= -tol_) & (abs(ov_) > tol_ | abs(bv_) > tol_);
        canonj   = nonzeromodesic( is_canon(nonzeromodesic));
        derivj   = nonzeromodesic(~is_canon(nonzeromodesic));
        runorder = [canonj, derivj];   % canonical first, so images can read them back
    else
        runorder = flip(nonzeromodesic);
    end

    for j = runorder
        if sym_inflow
            do_clst = is_canon(j);
        else
            do_clst = j >= round(Mode.nf/2);
            j2      = round(Mode.nf/2)-(j-round(Mode.nf/2));
        end

        if do_clst   % Physical / canonical mode -> CLST or linear-PSE profile
            if is_pse
                % Linear-PSE (non-parallel) seed: march DeHNSSo's own linear PSE from
                % upstream to R0 and take the developed shape (unit-max|u|) at R0,
                % carrying the mode's streamwise history instead of a local eigenfunction.
                % Compressible density-form (pse) or incompressible u,v,w,p form (pse_inc).
                if is_compr
                    [prof, al_pse] = solver_CPSE(ov_(j), bv_(j), StabGrid, Opt);
                    StabRes.rho(j,:,1)=prof.rho; StabRes.T(j,:,1)=prof.T;
                else
                    [prof, al_pse] = solver_IPSE(ov_(j), bv_(j), StabGrid, Opt);
                end
                StabRes.u(j,:,1)=prof.u; StabRes.v(j,:,1)=prof.v; StabRes.w(j,:,1)=prof.w; StabRes.p(j,:,1)=prof.p;
                StabRes.alpha(j,1)=al_pse;
            else
                StabRes = clst_inflow(j, 1, StabGrid, StabRes, bc_top, Opt);
            end

            % CLST/pse both return unit-max|u'|; scale by A0 and rebuild phi.
            StabRes.u(j,:,1)   = A0s(j).*StabRes.u(j,:,1);
            StabRes.v(j,:,1)   = A0s(j).*StabRes.v(j,:,1);
            StabRes.w(j,:,1)   = A0s(j).*StabRes.w(j,:,1);
            StabRes.p(j,:,1)   = A0s(j).*StabRes.p(j,:,1);
            if is_compr
                StabRes.rho(j,:,1) = A0s(j).*StabRes.rho(j,:,1);
                StabRes.T(j,:,1)   = A0s(j).*StabRes.T(j,:,1);
                StabRes.phi(j,:,1) = [StabRes.u(j,:,1), StabRes.v(j,:,1), ...
                                      StabRes.w(j,:,1), StabRes.rho(j,:,1), StabRes.T(j,:,1)];
            else
                StabRes.phi(j,:,1) = [StabRes.u(j,:,1), StabRes.v(j,:,1), ...
                                      StabRes.w(j,:,1), StabRes.p(j,:,1)];
            end

            switch lower(Stab.phaseRef)
                case 'pwall'
                    theta0 = angle(StabRes.p(j,end,1));
                case 'umax'
                    [~,iy] = max(abs(StabRes.u(j,:,1)));
                    theta0 = angle(StabRes.u(j,iy,1));
                otherwise
                    error('Unknown phaseRef option');
            end
            StabRes.phi(j,:,1) = StabRes.phi(j,:,1) * exp(-1i*theta0);

        elseif sym_inflow   % Symmetry image of a canonical mode (2D, Wbar=0)
            jc = find(abs(abs(ov_)-abs(ov_(j))) < tol_ & abs(abs(bv_)-abs(bv_(j))) < tol_ ...
                      & (ov_ >= -tol_) & (bv_ >= -tol_) & (abs(ov_) > tol_ | abs(bv_) > tol_), 1);
            sc      = Stab.A0(j)/Stab.A0(jc);   % rescale to this mode's amplitude
            do_conj = ov_(j) < -tol_;           % reality conjugation (omega sign flip)
            % Spanwise mirror flips w'; reality conjugation already flips beta's sign,
            % so mirror only when the target beta sign disagrees with it (XOR do_conj).
            % Otherwise the (-1,1) partner gets +conj(w) and the pair is non-physical.
            flip_w  = (bv_(j) < -tol_) ~= do_conj;
            if do_conj
                u=conj(StabRes.u(jc,:,1));   v=conj(StabRes.v(jc,:,1));   w=conj(StabRes.w(jc,:,1));
                p=conj(StabRes.p(jc,:,1));
                if is_compr; r=conj(StabRes.rho(jc,:,1)); T=conj(StabRes.T(jc,:,1)); end
                StabRes.alpha(j,1) = -conj(StabRes.alpha(jc,1));
            else
                u=StabRes.u(jc,:,1);   v=StabRes.v(jc,:,1);   w=StabRes.w(jc,:,1);
                p=StabRes.p(jc,:,1);
                if is_compr; r=StabRes.rho(jc,:,1); T=StabRes.T(jc,:,1); end
                StabRes.alpha(j,1) =  StabRes.alpha(jc,1);
            end
            if flip_w, w = -w; end
            StabRes.u(j,:,1)=sc*u; StabRes.v(j,:,1)=sc*v; StabRes.w(j,:,1)=sc*w; StabRes.p(j,:,1)=sc*p;
            if is_compr
                StabRes.rho(j,:,1)=sc*r; StabRes.T(j,:,1)=sc*T;
                StabRes.phi(j,:,1) = [StabRes.u(j,:,1), StabRes.v(j,:,1), ...
                                      StabRes.w(j,:,1), StabRes.rho(j,:,1), StabRes.T(j,:,1)];
            else
                StabRes.phi(j,:,1) = [StabRes.u(j,:,1), StabRes.v(j,:,1), ...
                                      StabRes.w(j,:,1), StabRes.p(j,:,1)];
            end

        else   % Conjugate mode (1D index-mirror path)
            StabRes.u(j,:,1)   = conj(StabRes.u(j2,:,1));
            StabRes.v(j,:,1)   = conj(StabRes.v(j2,:,1));
            StabRes.w(j,:,1)   = conj(StabRes.w(j2,:,1));
            StabRes.p(j,:,1)   = conj(StabRes.p(j2,:,1));
            if is_compr
                StabRes.rho(j,:,1) = conj(StabRes.rho(j2,:,1));
                StabRes.T(j,:,1)   = conj(StabRes.T(j2,:,1));
            end
            StabRes.alpha(j,1) = -conj(StabRes.alpha(j2,1));
            StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
        end
    end
    if is_pse
        fprintf('Initialising with linear-PSE (non-parallel) solution at inflow.\n')
    else
        fprintf('Initialising with CLST solution at inflow.\n')
    end

elseif Stab.IC=="ZERO"
    % No inflow disturbance: the field is generated entirely by wall/freestream
    % forcing (Stab.bcw / Stab.bct). StabRes.u/v/w/p (+ rho/T) keep their zero
    % initialisation at the inlet, and nonzeromodesic stays empty, so the active
    % modes come from nonzeromodeswall only (e.g. a wall actuator, roughness, or
    % the resonator impedance BC with no incoming instability wave).
    if isempty(nonzeromodeswall)
        warning('inflow:zeroNoForcing', ...
            ['Stab.IC=''ZERO'' with no wall/freestream forcing (Stab.bcw/bct all ' ...
             'zero): the disturbance field will be identically zero.']);
    end
    fprintf('Initialising with ZERO inflow (disturbance from wall/freestream forcing only).\n')

else
    error('Unknown inflow condition Stab.IC = ''%s''. Valid options: ''ILST'', ''CLST'', ''PSE'', ''LOAD'', ''ZERO''.', Stab.IC);
end

% PSE inflow override (Opt.pse_inflow): impose a developed multi-mode NPSE state
% at the inlet (ALL modes, not just the linear pair) so generated modes don't
% re-establish inside the domain. File fields: m,n, eta, u,v,w,rho,T,p (physical u_hat).
if isfield(Opt,'pse_inflow') && ~isempty(Opt.pse_inflow)
    P = load(Opt.pse_inflow);
    ovv = real(StabRes.omegavec(:)); bvv = real(StabRes.betavec(:)); tolp = 1e-9;
    o0 = min(abs(ovv(abs(ovv)>tolp))); b0 = min(abs(bvv(abs(bvv)>tolp)));
    if isempty(o0)||o0==0; o0=1; end; if isempty(b0)||b0==0; b0=1; end
    yq = StabGrid.etaun(:).'; loaded = [];
    for j = 1:Mode.nf
        mj = round(ovv(j)/o0); nj = round(bvv(j)/b0);
        k = find(P.m==mj & P.n==nj, 1); cc = false;
        if isempty(k); k = find(P.m==-mj & P.n==-nj, 1); cc = true; end   % reality conjugate
        if isempty(k); continue; end
        gv = @(F) interp1(P.eta(:).', F(k,:), yq, 'pchip', 0);
        u=gv(P.u); v=gv(P.v); w=gv(P.w); r=gv(P.rho); T=gv(P.T); p=gv(P.p);
        if cc; u=conj(u); v=conj(v); w=conj(w); r=conj(r); T=conj(T); p=conj(p); end
        StabRes.u(j,:,1)=u; StabRes.v(j,:,1)=v; StabRes.w(j,:,1)=w;
        StabRes.rho(j,:,1)=r; StabRes.T(j,:,1)=T; StabRes.p(j,:,1)=p;
        StabRes.phi(j,:,1)=[u,v,w,r,T]; loaded=[loaded j];
    end
    nonzeromodesic = unique([nonzeromodesic, loaded]);
    fprintf('PSE inflow: loaded %d modes from %s\n', numel(loaded), Opt.pse_inflow);
end

% nonzero modes are those initiated at inflow + those forced at wall
nonzeromodes = unique([nonzeromodesic, nonzeromodeswall]);
Mode.RunJ = nonzeromodes;            % Used in loops to assess only active modes
Mode.Aactive(Mode.RunJ) = 1;      % Flag these modes as active

% Define phiIC
StabRes.phiIC = squeeze(StabRes.phi(:,:,1));

end

%% ========================================================================
%  Per-mode inflow seeding: solve the local eigenvalue problem (solver_ILST /
%  solver_CLST), select the physical mode (EVfilter), normalise by max|u|, and
%  store the shape in StabRes. The eigen-solvers themselves live in init/lst/.
%% ========================================================================

function [StabRes] = ilst_inflow(j,i,StabGrid,StabRes,bc_top,Opt)
%ILST_INFLOW  Incompressible LST seed for mode j (Orr-Sommerfeld/Squire).

Re = StabGrid.Re;
U = StabGrid.U(:,i);
W = StabGrid.W(:,i);
dyU = StabGrid.dyU(:,i);
dyW = StabGrid.dyW(:,i);

omega = StabRes.omegavec(j);
beta = StabRes.betavec(j);
ny = StabGrid.ny;

if ~exist('bc_top','var'); bc_top = 'H_DR'; end % seed mode decays: freestream BC immaterial
if ~exist('Opt','var');   Opt.plot = 'on'; end
if ~isfield(Opt,'verbose'); Opt.verbose = 'on'; end

assert(isfield(StabGrid,'H') && isfield(StabGrid,'y_i'), 'ILST:gridMeta', ...
    'StabGrid needs H and y_i (regenerate the grid via gridgen).');

if strcmpi(Opt.verbose,'on')
    fprintf(' Calculating IC for nonzero modes with ILST: ')
end
% solver_ILST rebuilds the StabGrid.y grid (same H, y_i, ny) and returns the
% eigenfunctions on it, so no re-interpolation is needed on Malik grids.
[EV,u,v,w,p,yOS] = solver_ILST(Re,U,W,dyU,dyW,StabGrid.y,omega,beta,ny,StabGrid.H,StabGrid.y_i);

k = find(U<=0.99*max(U),1,'first');
d99 = StabGrid.y(k);
[alpha, index] = EVfilter(EV, v, yOS, omega, beta, d99);

if strcmpi(Opt.verbose,'on')
    fprintf('alpha = %.6f %+.6fi\n', real(alpha), imag(alpha));
end

StabRes.alpha(j,i) = alpha;

% Select the physical mode and put it on StabGrid.y (identity for Malik grids
% where yOS == StabGrid.y; a genuine interpolation for FD grids).
u = interp1(yOS, u(:,index), StabGrid.y, 'spline');
v = interp1(yOS, v(:,index), StabGrid.y, 'spline');
w = interp1(yOS, w(:,index), StabGrid.y, 'spline');
p = interp1(yOS, p(:,index), StabGrid.y, 'spline');

% Amplitude = max|u| after fine spline interpolation
y_mdInt     = linspace(StabGrid.y(1),StabGrid.y(end),4000);
u1int       = interp1(StabGrid.y,abs(u),y_mdInt,'spline');
[ampltd, ~] = max(abs(u1int));

StabRes.u(j,:,i)  = u/ampltd;
StabRes.v(j,:,i)  = v/ampltd;
StabRes.w(j,:,i)  = w/ampltd;
StabRes.p(j,:,i)  = p/ampltd;
StabRes.phi(j,:,i) = [ StabRes.u(j,:,i).'; StabRes.v(j,:,i).'; ...
                       StabRes.w(j,:,i).'; StabRes.p(j,:,i).'];

%% Visual check for mode selection
if ~strcmpi(Opt.plot, 'on'); return; end

figure(1); clf; set(gcf, 'Position', [200 200 1200 350]);
tiledlayout(1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
plot(U, StabGrid.y, 'b-', 'LineWidth', 1.5); hold on;
plot(W, StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(dyU, StabGrid.y, 'b--', 'LineWidth', 1.0);
plot(dyW, StabGrid.y, 'r--', 'LineWidth', 1.0);
xlabel('Velocity / derivative', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('U', 'W', 'dU/d\eta', 'dW/d\eta', 'Location', 'best', 'FontSize', 9)
title('Inflow base flow', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(real(EV), imag(EV), 'k.', 'MarkerSize', 12); hold on;
plot(real(alpha), imag(alpha), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
alpha_r_center = omega;  % entropy branch at c=1, so alpha_r = omega
xlim([alpha_r_center - 0.1, alpha_r_center + 0.1])
ylim([-0.1, 0.02])
xlabel('$\alpha_r$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\alpha_i$', 'Interpreter', 'latex', 'FontSize', 12)
legend({'Eigenvalues', 'Selected'}, 'Location', 'best', 'FontSize', 9)
title('Eigenvalue spectrum', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(abs(squeeze(StabRes.u(j,:,1))), StabGrid.y, 'k-', 'LineWidth', 1.5);
xlabel('$|\hat{u}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
title('Eigenfunction', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

sgtitle(sprintf('Inflow ILST: $\\alpha = %.4f %+.4fi$', real(alpha), imag(alpha)), ...
    'Interpreter', 'latex', 'FontSize', 14);
drawnow;

end

function [StabRes] = clst_inflow(j, i, StabGrid, StabRes, bc_top, Opt)
%CLST_INFLOW  Compressible LST seed for mode j, variables (rho,u,v,w,T).

ny    = StabGrid.ny;
Ma    = StabGrid.Ma;
ga    = StabGrid.gamma;
gM2   = ga * Ma^2;
N5    = 5 * ny;

omega = StabRes.omegavec(j);
beta  = StabRes.betavec(j);

if ~exist('bc_top','var'); bc_top = 'H_DR'; end
if ~exist('Opt','var');    Opt = struct('plot','on'); end
if ~isfield(Opt,'verbose'); Opt.verbose = 'on'; end

if strcmpi(Opt.verbose,'on')
    fprintf(' Calculating IC for nonzero modes with CLST: ')
end

U  = StabGrid.U(:,i);          % for d99 / plotting
Rb = StabGrid.rho(:,i);        % for p reconstruction (EOS)
Tb = StabGrid.T(:,i);

% Assemble and solve the compressible LST eigenvalue problem
drop_a2 = isfield(Opt,'drop_alpha2') && Opt.drop_alpha2;
[Eval, Evec] = solver_CLST(StabGrid, i, omega, beta, bc_top, drop_a2);

% Select the physical mode
if isfield(Opt, 'alpha_shift') && ~isempty(Opt.alpha_shift)
    % Continuation: nearest eigenvalue to the previous-station alpha, within a
    % physical band (avoids spurious / continuous-spectrum modes -> de-spikes sweep).
    band = real(Eval) > 0.05 & real(Eval) < 0.5 & abs(imag(Eval)) < 0.05;
    d = abs(Eval - Opt.alpha_shift);  d(~band) = inf;
    [~, index] = min(d);
    alpha = Eval(index);
else
    % Full spectrum + EVfilter (station 1: mode identification)
    v_evec = Evec(3:5:N5, :);
    kd99 = find(U <= 0.99*max(U), 1, 'first');
    d99_val = StabGrid.y(kd99);
    % Optional Mack-mode selector (Ma>1): Opt.mode_select = 'auto'|'slow'|'fast'
    mode_select = 'auto';
    if isfield(Opt,'mode_select') && ~isempty(Opt.mode_select)
        mode_select = Opt.mode_select;
    end
    [alpha, index] = EVfilter(Eval, v_evec, StabGrid.y, omega, beta, d99_val, Ma, mode_select);
end

StabRes.alpha(j,i) = alpha;
if strcmpi(Opt.verbose,'on')
    fprintf('alpha = %.6f %+.6fi\n', real(alpha), imag(alpha));
end

% Extract and normalise eigenfunction. CLST solves in (rho,u,v,w,T); convert to
% (u,v,w,p,T) via EOS: p' = (rho'*T_bar + rho_bar*T') / (gamma*Ma^2).
q       = Evec(1:N5, index);
rho_hat = q(1:5:N5);
u_hat   = q(2:5:N5);
v_hat   = q(3:5:N5);
w_hat   = q(4:5:N5);
T_hat   = q(5:5:N5);
p_hat   = (rho_hat .* Tb + Rb .* T_hat) / gM2;

y_int  = linspace(StabGrid.y(1), StabGrid.y(end), 4000);
u_int  = interp1(StabGrid.y, abs(u_hat), y_int, 'spline');
ampltd = max(abs(u_int));

StabRes.u(j,:,i)   = u_hat.'   / ampltd;
StabRes.v(j,:,i)   = v_hat.'   / ampltd;
StabRes.w(j,:,i)   = w_hat.'   / ampltd;
StabRes.p(j,:,i)   = p_hat.'   / ampltd;
StabRes.rho(j,:,i) = rho_hat.' / ampltd;

% Compressible LHS is density form: slot 4 of phi is rho'.
if isfield(StabRes,'T')
    StabRes.T(j,:,i) = T_hat.' / ampltd;
    StabRes.phi(j,:,i) = [StabRes.u(j,:,i),    StabRes.v(j,:,i), ...
                           StabRes.w(j,:,i),    StabRes.rho(j,:,i), ...
                           StabRes.T(j,:,i)];
else
    StabRes.phi(j,:,i) = [StabRes.u(j,:,i),  StabRes.v(j,:,i), ...
                           StabRes.w(j,:,i),  StabRes.rho(j,:,i)];
end

%% Plotting
if ~strcmpi(Opt.plot, 'on'); return; end

figure(1); clf; set(gcf, 'Position', [200 200 1400 350]);
tiledlayout(1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
% delta99: y where U first drops below 0.99*U_max (works ascending OR descending y)
kd99_plot    = find(U >= 0.99*max(U), 1, 'last');
if isempty(kd99_plot); kd99_plot = numel(U); end
delta99_plot = StabGrid.y(kd99_plot);
eta_max_plot = min(2 * delta99_plot, max(StabGrid.y));
plot(U, StabGrid.y, 'b-', 'LineWidth', 1.5); hold on;
plot(Tb, StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(Rb, StabGrid.y, 'g-', 'LineWidth', 1.5);
ylim([0, eta_max_plot]);
xlabel('Base flow', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$U$', '$T$', '$\rho$', 'Location', 'best', 'Interpreter', 'latex')
title('Inflow base flow', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
phys = ~isinf(Eval) & ~isnan(Eval) & abs(Eval) < 10;
plot(real(Eval(phys)), imag(Eval(phys)), 'k.', 'MarkerSize', 12); hold on;
plot(real(alpha), imag(alpha), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
xlim([real(alpha)-0.2, real(alpha)+0.2]);
ylim([-0.1, 0.05]);
xlabel('$\alpha_r$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\alpha_i$', 'Interpreter', 'latex', 'FontSize', 12)
legend({'Eigenvalues', 'Selected'}, 'Location', 'best', 'FontSize', 9)
title('Eigenvalue spectrum', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(abs(squeeze(StabRes.u(j,:,i))), StabGrid.y, 'k-', 'LineWidth', 1.5); hold on;
plot(abs(squeeze(StabRes.v(j,:,i))), StabGrid.y, 'b--', 'LineWidth', 1.2);
xlabel('$|\hat{q}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
legend('$|\hat{u}|$', '$|\hat{v}|$', 'Interpreter', 'latex', 'Location', 'best')
title('Velocity eigenfunctions', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

nexttile
plot(abs(squeeze(StabRes.p(j,:,i))), StabGrid.y, 'g-', 'LineWidth', 1.5); hold on;
plot(abs(T_hat)/ampltd,               StabGrid.y, 'r-', 'LineWidth', 1.5);
plot(abs(rho_hat)/ampltd,             StabGrid.y, 'm-', 'LineWidth', 1.5);
legend({'$|\hat{p}|$', '$|\hat{T}|$', '$|\hat{\rho}|$'}, ...
       'Interpreter', 'latex', 'Location', 'best')
xlabel('$|\hat{q}|$', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('$\eta$', 'Interpreter', 'latex', 'FontSize', 12)
title('Thermodynamic eigenfunctions', 'Interpreter', 'latex', 'FontSize', 13)
grid on; box on;

end

