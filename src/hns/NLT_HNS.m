function [ f ] = NLT_HNS(Grid,RunJ,StabRes,D1,HB)
%NLT_HNS  Mode-to-mode nonlinear forcing for the HNS system.
%
%   f = NLT_HNS(Grid, RunJ, StabRes, D1, HB)
%
%   Evaluates the nonlinear forcing that couples the harmonic modes of the
%   harmonic linearised Navier-Stokes (HNS) solver. For each active mode j the
%   harmonic-balancing registry HB selects the interacting mode pairs (and, in
%   the compressible case, triples); their products build the momentum,
%   continuity and energy forcing. Incompressible runs (nvar=4) force momentum
%   only. Compressible runs (nvar=5) add the density-weighted convection, the
%   nonlinear EOS pressure, the rho'-coupling, variable-viscosity and
%   pressure-dilatation groups, consistent with LHS_compressible (c_v form).
%
%   Inputs:
%     Grid     — grid, metrics and base-flow fields
%     RunJ     — indices of the active modes
%     StabRes  — current mode fields (u, v, w, p, rho, T, q)
%     D1       — wall-normal first-derivative operator
%     HB       — harmonic-balancing interaction registry (see Hbalancing)
%   Outputs:
%     f        — [nf x nvar*ny*nx] nonlinear forcing per mode

%% Find dimensions
[nf,ny,nx] = size(StabRes.u);
iu = 1i;

if isfield(StabRes, 'q') && ~isempty(StabRes.q)
    nvar = round(size(StabRes.q, 2) / (ny*nx));   % 4 incompressible | 5 compressible
elseif isfield(StabRes, 'T') && ~isempty(StabRes.T)
    nvar = 5;
else
    nvar = 4;
end
f    = zeros(nf, nvar*ny*nx);

% Metric coefficients
etay = Grid.etay;
etax = Grid.etax;
xiy  = Grid.xiy;
xix  = Grid.xix;

% ξ-derivatives via sparse non-uniform FD matrix (built once in init.m).
% On a uniform grid this reduces to classical 4th-order coefficients; on
% xrefined / curvilinear grids it correctly accounts for variable Δξ.
D_xi = Grid.D_xi;

%% Calculate derivative fields (batched over modes)

% xi-derivatives: apply the sparse D_xi matrix along rows (streamwise)
dudxi = zeros(nf, ny, nx);
dvdxi = zeros(nf, ny, nx);
dwdxi = zeros(nf, ny, nx);
for j = RunJ
    dudxi(j,:,:) = squeeze(StabRes.u(j,:,:)) * D_xi.';
    dvdxi(j,:,:) = squeeze(StabRes.v(j,:,:)) * D_xi.';
    dwdxi(j,:,:) = squeeze(StabRes.w(j,:,:)) * D_xi.';
end

% eta-derivatives: apply D1 to all stations at once per mode
dudeta = zeros(nf, ny, nx);
dvdeta = zeros(nf, ny, nx);
dwdeta = zeros(nf, ny, nx);
for j = RunJ
    dudeta(j,:,:) = D1 * squeeze(StabRes.u(j,:,:));
    dvdeta(j,:,:) = D1 * squeeze(StabRes.v(j,:,:));
    dwdeta(j,:,:) = D1 * squeeze(StabRes.w(j,:,:));
end

% Compressible: density / temperature / pressure derivatives for the continuity
% (-div(rho'u')) and energy (-rho_bar*(u'.grad)T' + Ec*(u'.grad)p') NLT.
if nvar == 5
    drdxi  = zeros(nf, ny, nx);  drdeta = zeros(nf, ny, nx);
    dTdxi  = zeros(nf, ny, nx);  dTdeta = zeros(nf, ny, nx);
    dpdxi  = zeros(nf, ny, nx);  dpdeta = zeros(nf, ny, nx);
    for j = RunJ
        drdxi(j,:,:)  = squeeze(StabRes.rho(j,:,:)) * D_xi.';
        drdeta(j,:,:) = D1 * squeeze(StabRes.rho(j,:,:));
        dTdxi(j,:,:)  = squeeze(StabRes.T(j,:,:)) * D_xi.';
        dTdeta(j,:,:) = D1 * squeeze(StabRes.T(j,:,:));
        dpdxi(j,:,:)  = squeeze(StabRes.p(j,:,:)) * D_xi.';
        dpdeta(j,:,:) = D1 * squeeze(StabRes.p(j,:,:));
    end
    Ec = (Grid.gamma - 1) * Grid.Ma^2;        % Eckert number = (gamma-1) Ma^2
    % base-flow fields/gradients for the rho'-coupling group -rho'(d_t + Ubar.grad + u'.grad Ubar)
    Ubar = Grid.U; Vbar = Grid.V; Wbar = Grid.W;
    bdxU = Grid.dxU; bdyU = Grid.dyU; bdxV = Grid.dxV; bdyV = Grid.dyV;
    bdxW = Grid.dxW; bdyW = Grid.dyW; bdyT = Grid.dyT; bdxT = Grid.dxT;
    % base spanwise (d_z) gradients: 0 for a quasi-2D base, nonzero for a 3D/swept base
    bdzU = zeros(ny,nx); bdzV = zeros(ny,nx); bdzW = zeros(ny,nx); bdzT = zeros(ny,nx);
    if isfield(Grid,'dzU')&&~isempty(Grid.dzU), bdzU = Grid.dzU; end
    if isfield(Grid,'dzV')&&~isempty(Grid.dzV), bdzV = Grid.dzV; end
    if isfield(Grid,'dzW')&&~isempty(Grid.dzW), bdzW = Grid.dzW; end
    if isfield(Grid,'dzT')&&~isempty(Grid.dzT), bdzT = Grid.dzT; end
    % cubic (3-wave) interaction lookup: spectral (m,n) index -> mode number
    ov=real(StabRes.omegavec(:)).'; bv=real(StabRes.betavec(:)).';
    om0c=min(abs(ov(abs(ov)>1e-30))); if isempty(om0c), om0c=1; end
    bt0c=min(abs(bv(abs(bv)>1e-30))); if isempty(bt0c), bt0c=1; end
    mmv=round(ov/om0c); nnv=round(bv/bt0c);
    mlo=min(mmv); mhi=max(mmv); nlo=min(nnv); nhi=max(nnv);
    Lkp=zeros(mhi-mlo+1, nhi-nlo+1);
    for kk=1:nf, Lkp(mmv(kk)-mlo+1, nnv(kk)-nlo+1)=kk; end
end

%% Calculate forcing (local ny x nx accumulation per mode, avoids 3D indexing)

for j = RunJ(RunJ >= round(nf/2))

    % Find interacting mode pairs from Harmonic Balancing
    [gharray, jkarray] = find(HB(:,:,j));

    % Keep only active mode pairs
    active = ismember(gharray, RunJ) & ismember(jkarray, RunJ);
    gharray = gharray(active);
    jkarray = jkarray(active);

    % Local accumulation arrays — avoids squeeze into 3D array each interaction
    xm_j = zeros(ny, nx);
    ym_j = zeros(ny, nx);
    zm_j = zeros(ny, nx);
    cm_j = zeros(ny, nx);          % continuity (rho) NLT: -div(rho' u')
    em_j = zeros(ny, nx);          % energy (T) NLT: -rho_bar*(u'.grad)T' + Ec*(u'.grad)p'
    xmr_j = zeros(ny, nx);         % momentum rho'-coupling (not rho_bar-weighted)
    ymr_j = zeros(ny, nx);
    zmr_j = zeros(ny, nx);

    for j2 = 1:length(gharray)
        gh = gharray(j2);
        jk = jkarray(j2);
        beta = StabRes.betavec(jk);

        % Extract fields for this mode pair (squeeze once)
        u_gh = squeeze(StabRes.u(gh,:,:));
        v_gh = squeeze(StabRes.v(gh,:,:));
        w_gh = squeeze(StabRes.w(gh,:,:));

        dudxi_jk  = squeeze(dudxi(jk,:,:));
        dvdxi_jk  = squeeze(dvdxi(jk,:,:));
        dwdxi_jk  = squeeze(dwdxi(jk,:,:));
        dudeta_jk = squeeze(dudeta(jk,:,:));
        dvdeta_jk = squeeze(dvdeta(jk,:,:));
        dwdeta_jk = squeeze(dwdeta(jk,:,:));
        u_jk = squeeze(StabRes.u(jk,:,:));
        v_jk = squeeze(StabRes.v(jk,:,:));
        w_jk = squeeze(StabRes.w(jk,:,:));

        % Accumulate momentum forcing into local arrays
        xm_j = xm_j ...
            - u_gh.*dudxi_jk.*xix - u_gh.*dudeta_jk.*etax ...
            - v_gh.*dudxi_jk.*xiy - v_gh.*dudeta_jk.*etay ...
            - iu*beta .* w_gh.*u_jk;

        ym_j = ym_j ...
            - u_gh.*dvdxi_jk.*xix - u_gh.*dvdeta_jk.*etax ...
            - v_gh.*dvdxi_jk.*xiy - v_gh.*dvdeta_jk.*etay ...
            - iu*beta .* w_gh.*v_jk;

        zm_j = zm_j ...
            - u_gh.*dwdxi_jk.*xix - u_gh.*dwdeta_jk.*etax ...
            - v_gh.*dwdxi_jk.*xiy - v_gh.*dwdeta_jk.*etay ...
            - iu*beta .* w_gh.*w_jk;

        % Continuity NLT: -div(rho' u') = -[d/dx(rho_gh u_jk) + d/dy(rho_gh v_jk)
        % + d/dz(rho_gh w_jk)], product rule with each factor carrying its wavenumber.
        if nvar == 5
            rho_gh    = squeeze(StabRes.rho(gh,:,:));
            beta_gh   = StabRes.betavec(gh);
            drdxi_gh  = squeeze(drdxi(gh,:,:));
            drdeta_gh = squeeze(drdeta(gh,:,:));
            drho_dx = drdxi_gh.*xix + drdeta_gh.*etax;
            drho_dy = drdxi_gh.*xiy + drdeta_gh.*etay;
            du_jk_dx = dudxi_jk.*xix + dudeta_jk.*etax;
            dv_jk_dy = dvdxi_jk.*xiy + dvdeta_jk.*etay;
            cm_j = cm_j ...
                - ( drho_dx.*u_jk + rho_gh.*du_jk_dx ) ...
                - ( drho_dy.*v_jk + rho_gh.*dv_jk_dy ) ...
                - ( iu*beta_gh.*rho_gh.*w_jk + rho_gh.*(iu*beta.*w_jk) );

            % Energy (c_v form, consistent with LHS_compressible): convection
            % -rho_bar*(u'.grad)T'  + pressure-DILATATION -(gamma-1)(rho T div u)'
            % [added below where div u'_jk is available]. Was cp/pressure-work +Ec*(u'.grad)p'.
            T_jk = squeeze(StabRes.T(jk,:,:));
            dT_jk_dx = squeeze(dTdxi(jk,:,:)).*xix + squeeze(dTdeta(jk,:,:)).*etax;
            dT_jk_dy = squeeze(dTdxi(jk,:,:)).*xiy + squeeze(dTdeta(jk,:,:)).*etay;
            convT = u_gh.*dT_jk_dx + v_gh.*dT_jk_dy + iu*beta.*w_gh.*T_jk;
            em_j = em_j - Grid.rho.*convT;

            % Nonlinear EOS pressure -grad(rho'T')/(gamma Ma^2) (quadratic): the
            % reconstructed p' is linear in (rho',T'), so the rho'T' piece is
            % absent from the density-form momentum gradient. From CNPSE nlt_forcing.
            % Goes in xmr_j (NOT rho_bar-weighted; xm_j is multiplied by rho_b below).
            ikM2 = 1/(Grid.gamma*Grid.Ma^2);
            xmr_j = xmr_j - ikM2*( drho_dx.*T_jk + rho_gh.*dT_jk_dx );
            ymr_j = ymr_j - ikM2*( drho_dy.*T_jk + rho_gh.*dT_jk_dy );
            zmr_j = zmr_j - ikM2*( iu*beta_gh.*rho_gh.*T_jk + rho_gh.*(iu*beta.*T_jk) );

            % (c_v form: the cp EOS pressure-work +Ec*D(rho'T')/Dt/gM2 is dropped;
            %  the pressure term is the dilatation added below.)

            % rho'-coupling group: -rho'( d_t q' + Ubar.grad q' + u'.grad Qbar ),
            % rho' = rho_gh, d_t -> -i*omega_jk, base d_z = 0.
            omega_jk = StabRes.omegavec(jk);
            du_jk_dy = squeeze(dudxi(jk,:,:)).*xiy + squeeze(dudeta(jk,:,:)).*etay;
            dv_jk_dx = squeeze(dvdxi(jk,:,:)).*xix + squeeze(dvdeta(jk,:,:)).*etax;
            dw_jk_dx = squeeze(dwdxi(jk,:,:)).*xix + squeeze(dwdeta(jk,:,:)).*etax;
            dw_jk_dy = squeeze(dwdxi(jk,:,:)).*xiy + squeeze(dwdeta(jk,:,:)).*etay;
            xmr_j = xmr_j - rho_gh.*( -iu*omega_jk.*u_jk + Ubar.*du_jk_dx + Vbar.*du_jk_dy ...
                          + Wbar.*(iu*beta.*u_jk) + u_jk.*bdxU + v_jk.*bdyU + w_jk.*bdzU );
            ymr_j = ymr_j - rho_gh.*( -iu*omega_jk.*v_jk + Ubar.*dv_jk_dx + Vbar.*dv_jk_dy ...
                          + Wbar.*(iu*beta.*v_jk) + u_jk.*bdxV + v_jk.*bdyV + w_jk.*bdzV );
            zmr_j = zmr_j - rho_gh.*( -iu*omega_jk.*w_jk + Ubar.*dw_jk_dx + Vbar.*dw_jk_dy ...
                          + Wbar.*(iu*beta.*w_jk) + u_jk.*bdxW + v_jk.*bdyW + w_jk.*bdzW );
            em_j  = em_j  - rho_gh.*( -iu*omega_jk.*T_jk + Ubar.*dT_jk_dx + Vbar.*dT_jk_dy ...
                          + Wbar.*(iu*beta.*T_jk) + u_jk.*bdxT + v_jk.*bdyT + w_jk.*bdzT );

            % --- variable-mu viscous group (flat-plate 2nd derivatives, as in CNPSE) ---
            % viscosity perturbation mu'(T)=muT*T_gh; lambda=-2/3 mu (Stokes), k=mu.
            muT_b = Grid.dmuT;  muTT_b = Grid.d2muT2;  mu_b = Grid.mu;
            % Stokes bulk viscosity lambda = -2/3 mu (consistent with LHS_compressible).
            % NB: CNPSE nlt_forcing uses lamm=mum (lambda=mu), inconsistent with its own
            % Stokes linear operator; a latent CNPSE bug, harmless there (bulk ~O(div u')~O(Ma^2)).
            lamT_b = -2/3*muT_b; lamTT_b = -2/3*muTT_b; lam_b = -2/3*mu_b;
            Re = Grid.Re;  RePr = Grid.Re*Grid.Pr;  T_gh = squeeze(StabRes.T(gh,:,:));
            % jk second physical derivatives via the full curvilinear chain rule, the
            % same expansion the base flow is differentiated with in compute_bf_derivatives:
            %   d_xx q = xix^2 q_xixi + 2 xix etax q_xieta + etax^2 q_etaeta + xixx q_xi + etaxx q_eta
            % and analogously for d_yy, d_xy (first AND second inverse metrics). The
            % computational parts use D2 (eta-eta), Grid.D_xi2 (xi-xi, 4th-order) and
            % D1 (eta of q_xi). Reduces to D_xi2 / D2 / D1.D_xi on a Cartesian grid.
            D2 = Grid.D2;
            if isfield(Grid,'D_xi2') && ~isempty(Grid.D_xi2), D_xi2 = Grid.D_xi2; else, D_xi2 = D_xi*D_xi; end
            mets = struct('xix',xix,'etax',etax,'xiy',xiy,'etay',etay, ...
                'xixx',Grid.xixx,'etaxx',Grid.etaxx,'xiyy',Grid.xiyy,'etayy',Grid.etayy, ...
                'xixy',Grid.xixy,'etaxy',Grid.etaxy);
            [uxx,uyy,uxy] = phys_d2(squeeze(StabRes.u(jk,:,:)), squeeze(dudxi(jk,:,:)), squeeze(dudeta(jk,:,:)), D_xi2,D1,D2,mets);
            [vxx,vyy,vxy] = phys_d2(squeeze(StabRes.v(jk,:,:)), squeeze(dvdxi(jk,:,:)), squeeze(dvdeta(jk,:,:)), D_xi2,D1,D2,mets);
            [wxx,wyy]     = phys_d2(squeeze(StabRes.w(jk,:,:)), squeeze(dwdxi(jk,:,:)), squeeze(dwdeta(jk,:,:)), D_xi2,D1,D2,mets);
            % gh first derivatives (for the strain-rate cross terms)
            dug_x = squeeze(dudxi(gh,:,:)).*xix + squeeze(dudeta(gh,:,:)).*etax;
            dug_y = squeeze(dudxi(gh,:,:)).*xiy + squeeze(dudeta(gh,:,:)).*etay;
            dvg_x = squeeze(dvdxi(gh,:,:)).*xix + squeeze(dvdeta(gh,:,:)).*etax;
            dvg_y = squeeze(dvdxi(gh,:,:)).*xiy + squeeze(dvdeta(gh,:,:)).*etay;
            dwg_y = squeeze(dwdxi(gh,:,:)).*xiy + squeeze(dwdeta(gh,:,:)).*etay;
            dTg_x = squeeze(dTdxi(gh,:,:)).*xix + squeeze(dTdeta(gh,:,:)).*etax;
            dTg_y = squeeze(dTdxi(gh,:,:)).*xiy + squeeze(dTdeta(gh,:,:)).*etay;
            divU = du_jk_dx + dv_jk_dy + iu*beta.*w_jk;     % div u'_jk
            % c_v pressure-dilatation -(gamma-1)(rho T div u)' [quadratic], matching
            % LHS_compressible; replaces the cp pressure-work removed above.
            em_j = em_j - (Grid.gamma-1)*( Grid.rho.*T_gh.*divU ...
                                         + rho_gh.*Grid.T.*divU ...
                                         + rho_gh.*T_jk.*(bdxU+bdyV+bdzW) );
            xmr_j = xmr_j ...
                + (1/Re)*lamT_b.*T_gh.*( uxx + vxy + iu*beta.*dw_jk_dx ) ...
                + (2/Re)*muT_b.*T_gh.*uxx ...
                + (1/Re)*muT_b.*T_gh.*( uyy + vxy - beta^2.*u_jk + iu*beta.*dw_jk_dx ) ...
                + (1/Re)*( lamTT_b.*bdxT.*T_gh + lamT_b.*dTg_x ).*divU ...
                + (2/Re)*( muTT_b.*bdxT.*T_gh + muT_b.*dTg_x ).*du_jk_dx ...
                + (1/Re)*( muTT_b.*bdyT.*T_gh + muT_b.*dTg_y ).*( du_jk_dy + dv_jk_dx ) ...
                + (1/Re)*( muTT_b.*bdzT.*T_gh + muT_b.*iu*beta_gh.*T_gh ).*( iu*beta.*u_jk + dw_jk_dx );
            ymr_j = ymr_j ...
                + (1/Re)*lamT_b.*T_gh.*( uxy + vyy + iu*beta.*dw_jk_dy ) ...
                + (1/Re)*muT_b.*T_gh.*( 2*vyy + uxy ) ...
                + (1/Re)*muT_b.*T_gh.*( vxx - beta^2.*v_jk + iu*beta.*dw_jk_dy ) ...
                + (1/Re)*( lamTT_b.*bdyT.*T_gh + lamT_b.*dTg_y ).*divU ...
                + (2/Re)*( muTT_b.*bdyT.*T_gh + muT_b.*dTg_y ).*dv_jk_dy ...
                + (1/Re)*( muTT_b.*bdxT.*T_gh + muT_b.*dTg_x ).*( du_jk_dy + dv_jk_dx ) ...
                + (1/Re)*( muTT_b.*bdzT.*T_gh + muT_b.*iu*beta_gh.*T_gh ).*( iu*beta.*v_jk + dw_jk_dy );
            zmr_j = zmr_j ...
                + (1/Re)*lamT_b.*T_gh.*( iu*beta.*du_jk_dx + iu*beta.*dv_jk_dy - beta^2.*w_jk ) ...
                + (1/Re)*muT_b.*T_gh.*( -2*beta^2.*w_jk + iu*beta.*du_jk_dx ) ...
                + (1/Re)*muT_b.*T_gh.*( wxx + iu*beta.*dv_jk_dy + wyy ) ...
                + (1/Re)*( lamTT_b.*bdzT.*T_gh + lamT_b.*iu*beta_gh.*T_gh ).*divU ...
                + (2/Re)*( muTT_b.*bdzT.*T_gh + muT_b.*iu*beta_gh.*T_gh ).*(iu*beta.*w_jk) ...
                + (1/Re)*( muTT_b.*bdxT.*T_gh + muT_b.*dTg_x ).*( iu*beta.*u_jk + dw_jk_dx ) ...
                + (1/Re)*( muTT_b.*bdyT.*T_gh + muT_b.*dTg_y ).*( iu*beta.*v_jk + dw_jk_dy );

            % --- energy viscous: heat conduction (k=mu) + viscous dissipation ---
            [Txx,Tyy] = phys_d2(squeeze(StabRes.T(jk,:,:)), squeeze(dTdxi(jk,:,:)), squeeze(dTdeta(jk,:,:)), D_xi2,D1,D2,mets);
            dwg_x = squeeze(dwdxi(gh,:,:)).*xix + squeeze(dwdeta(gh,:,:)).*etax;
            em_j = em_j ...
                + (Grid.gamma/RePr)*( muTT_b.*bdxT.*T_gh + muT_b.*dTg_x ).*dT_jk_dx ...
                + (Grid.gamma/RePr)*muT_b.*T_gh.*Txx ...
                + (Grid.gamma/RePr)*( muTT_b.*bdyT.*T_gh + muT_b.*dTg_y ).*dT_jk_dy ...
                + (Grid.gamma/RePr)*muT_b.*T_gh.*Tyy ...
                + (Grid.gamma/RePr)*( muTT_b.*bdzT.*T_gh + muT_b.*iu*beta_gh.*T_gh ).*( iu*beta.*T_jk ) ...
                + (Grid.gamma/RePr)*muT_b.*T_gh.*( -beta^2.*T_jk ) ...
                + (Grid.gamma*Ec/Re)*lam_b.*( dug_x.*du_jk_dx + dvg_y.*dv_jk_dy + (iu*beta_gh.*w_gh).*(iu*beta.*w_jk) ) ...
                + (Grid.gamma*Ec/Re)*lam_b.*( 2*dug_x.*(dv_jk_dy + iu*beta.*w_jk) + 2*dvg_y.*(iu*beta.*w_jk) ) ...
                + (Grid.gamma*Ec/Re)*lamT_b.*(2*(bdxU+bdyV+bdzW)).*T_gh.*divU ...
                + (2*Grid.gamma*Ec/Re)*mu_b.*( dug_x.*du_jk_dx + dvg_y.*dv_jk_dy + (iu*beta_gh.*w_gh).*(iu*beta.*w_jk) ) ...
                + (2*Grid.gamma*Ec/Re)*2*muT_b.*T_gh.*( bdxU.*du_jk_dx + bdyV.*dv_jk_dy + bdzW.*(iu*beta.*w_jk) ) ...
                + (Grid.gamma*Ec/Re)*mu_b.*( dug_y.*du_jk_dy + 2*dug_y.*dv_jk_dx + dvg_x.*dv_jk_dx ) ...
                + (Grid.gamma*Ec/Re)*2*muT_b.*T_gh.*( bdyU.*du_jk_dy + bdyU.*dv_jk_dx + bdxV.*du_jk_dy + bdxV.*dv_jk_dx ) ...
                + (Grid.gamma*Ec/Re)*mu_b.*( dwg_y.*dw_jk_dy + 2*dwg_y.*(iu*beta.*v_jk) + (iu*beta_gh.*v_gh).*(iu*beta.*v_jk) ) ...
                + (Grid.gamma*Ec/Re)*2*muT_b.*T_gh.*( bdyW.*dw_jk_dy + bdyW.*(iu*beta.*v_jk) ) ...
                + (Grid.gamma*Ec/Re)*mu_b.*( dwg_x.*dw_jk_dx + 2*dwg_x.*(iu*beta.*u_jk) + (iu*beta_gh.*u_gh).*(iu*beta.*u_jk) ) ...
                + (Grid.gamma*Ec/Re)*2*muT_b.*T_gh.*( bdxW.*(iu*beta.*w_jk + dw_jk_dx) + bdxW.*(iu*beta.*u_jk) );
        end
    end

    % --- cubic (3-wave) NLT: -rho'_a (u'_b.grad)q'_c (momentum/energy) + cubic
    % viscous dissipation. Triples (a,b,c) with spectral indices summing to mode j;
    % skip MFD-coupling triples unless all three are the MFD (matches CNPSE).
    if nvar == 5
        mj=mmv(j); nj=nnv(j);
        for a = RunJ
          for b = RunJ
            m3i=mj-mmv(a)-mmv(b); n3i=nj-nnv(a)-nnv(b);
            if m3i<mlo||m3i>mhi||n3i<nlo||n3i>nhi, continue; end
            c=Lkp(m3i-mlo+1, n3i-nlo+1);
            if c==0 || ~any(RunJ==c), continue; end
            % Full physics: cubic MFD<->wave triples retained (no skip).
            rho1=squeeze(StabRes.rho(a,:,:)); T1c=squeeze(StabRes.T(a,:,:));
            bb=StabRes.betavec(b); bc=StabRes.betavec(c);
            u2=squeeze(StabRes.u(b,:,:)); v2=squeeze(StabRes.v(b,:,:)); w2=squeeze(StabRes.w(b,:,:));
            u3=squeeze(StabRes.u(c,:,:)); v3=squeeze(StabRes.v(c,:,:)); w3=squeeze(StabRes.w(c,:,:)); T3=squeeze(StabRes.T(c,:,:));
            du3x=squeeze(dudxi(c,:,:)).*xix+squeeze(dudeta(c,:,:)).*etax; du3y=squeeze(dudxi(c,:,:)).*xiy+squeeze(dudeta(c,:,:)).*etay;
            dv3x=squeeze(dvdxi(c,:,:)).*xix+squeeze(dvdeta(c,:,:)).*etax; dv3y=squeeze(dvdxi(c,:,:)).*xiy+squeeze(dvdeta(c,:,:)).*etay;
            dw3x=squeeze(dwdxi(c,:,:)).*xix+squeeze(dwdeta(c,:,:)).*etax; dw3y=squeeze(dwdxi(c,:,:)).*xiy+squeeze(dwdeta(c,:,:)).*etay;
            dT3x=squeeze(dTdxi(c,:,:)).*xix+squeeze(dTdeta(c,:,:)).*etax; dT3y=squeeze(dTdxi(c,:,:)).*xiy+squeeze(dTdeta(c,:,:)).*etay;
            % cubic convection -rho'_a (u'_b.grad) q'_c
            xmr_j = xmr_j - rho1.*( u2.*du3x + v2.*du3y + iu*bc.*w2.*u3 );
            ymr_j = ymr_j - rho1.*( u2.*dv3x + v2.*dv3y + iu*bc.*w2.*v3 );
            zmr_j = zmr_j - rho1.*( u2.*dw3x + v2.*dw3y + iu*bc.*w2.*w3 );
            em_j  = em_j  - rho1.*( u2.*dT3x + v2.*dT3y + iu*bc.*w2.*T3 );
            % cubic dissipation: mu'(T'_a) * strain(b) * strain(c)
            du2x=squeeze(dudxi(b,:,:)).*xix+squeeze(dudeta(b,:,:)).*etax; du2y=squeeze(dudxi(b,:,:)).*xiy+squeeze(dudeta(b,:,:)).*etay;
            dv2x=squeeze(dvdxi(b,:,:)).*xix+squeeze(dvdeta(b,:,:)).*etax; dv2y=squeeze(dvdxi(b,:,:)).*xiy+squeeze(dvdeta(b,:,:)).*etay;
            dw2x=squeeze(dwdxi(b,:,:)).*xix+squeeze(dwdeta(b,:,:)).*etax; dw2y=squeeze(dwdxi(b,:,:)).*xiy+squeeze(dwdeta(b,:,:)).*etay;
            ibw2=iu*bb.*w2; ibv2=iu*bb.*v2; ibu2=iu*bb.*u2;
            ibw3=iu*bc.*w3; ibv3=iu*bc.*v3; ibu3=iu*bc.*u3;
            em_j = em_j ...
              + (Grid.gamma*Ec/Re)*lamT_b.*T1c.*( du2x.*du3x + dv2y.*dv3y + ibw2.*ibw3 + 2*du2x.*dv3y + 2*du2x.*ibw3 + 2*dv2y.*ibw3 ) ...
              + (2*Grid.gamma*Ec/Re)*muT_b.*T1c.*( du2x.*du3x + dv2y.*dv3y + ibw2.*ibw3 ) ...
              + (Grid.gamma*Ec/Re)*muT_b.*T1c.*( du2y.*du3y + 2*du2y.*dv3x + dv2x.*dv3x ) ...
              + (Grid.gamma*Ec/Re)*muT_b.*T1c.*( dw2y.*dw3y + 2*dw2y.*ibv3 + ibv2.*ibv3 ) ...
              + (Grid.gamma*Ec/Re)*muT_b.*T1c.*( dw2x.*dw3x + 2*dw2x.*ibu3 + ibu2.*ibu3 );
            % c_v cubic pressure-dilatation -(gamma-1) rho'_a T'_c (div u'_b):
            % replaces the cp cubic pressure-work +Ec*u'.grad(rho'T'/gM2).
            em_j = em_j - (Grid.gamma-1)*rho1.*T3.*( du2x + dv2y + iu*bb.*w2 );
          end
        end
    end

    % Compressible (nvar=5): density-weight the convective momentum forcing,
    % -rho_bar*(u'.grad)u'. At Ma->0 the base density rho_bar->1 so this reduces
    % to the incompressible forcing (recovery gate). Continuity (rho) row = cm_j
    % and energy (T) row = em_j are fully assembled above (quad + cubic).
    if nvar == 5
        rho_b = Grid.rho;                 % base density on the stab grid (ny x nx)
        mfmode = getenv('INCOMPF');       % momentum-forcing diagnostic
        mfsc = getenv('MFSCALE'); if isempty(mfsc), mfsc=1; else, mfsc=str2double(mfsc); end
        if isempty(mfmode)                % default: full compressible (mfsc=1); MFSCALE scales it
            xm_j = mfsc*(rho_b.*xm_j + xmr_j);  ym_j = mfsc*(rho_b.*ym_j + ymr_j);  zm_j = mfsc*(rho_b.*zm_j + zmr_j);
        elseif strcmp(mfmode,'rhob')      % rho_bar weight only (no rho'-coupling)
            xm_j = rho_b.*xm_j;          ym_j = rho_b.*ym_j;          zm_j = rho_b.*zm_j;
        elseif strcmp(mfmode,'xmr')       % rho'-coupling only (no rho_bar weight)
            xm_j = xm_j + xmr_j;         ym_j = ym_j + ymr_j;         zm_j = zm_j + zmr_j;
        end                               % else (INCOMPF=1): plain -(u'.grad)u' incompressible momentum forcing
        block = [xm_j(:,2:nx); ym_j(:,2:nx); zm_j(:,2:nx); cm_j(:,2:nx); em_j(:,2:nx)];
    else
        block = [xm_j(:,2:nx); ym_j(:,2:nx); zm_j(:,2:nx); zeros(ny, nx-1)];
    end
    f(j, nvar*ny+1:end) = block(:).';

end

end
