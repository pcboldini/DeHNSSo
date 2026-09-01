function F_real = newton_residual(q_real, ctx)
%NEWTON_RESIDUAL  Residual F(q) of the nonlinear HNS system, in real arithmetic.
%
%   F_real = newton_residual(q_real, ctx)
%
%   The state is real, q_real = [real(q_phys); imag(q_phys)], of length
%   2*nphys*bsize over the physical modes only; the conjugate modes are rebuilt
%   inside. Real arithmetic is what makes the Newton-Krylov linearisation well
%   posed: the NLT depends on both q_phys and conj(q_phys), so the complex
%   Jacobian is antilinear, and splitting into real and imaginary parts gives
%   GMRES the genuinely linear operator it needs.

%% Unpack state
% jc = nf+1-j pairs each physical mode with its conjugate slot; NLT_HNS needs the
% whole spectrum, so both halves are filled even though only phys is solved for.
nphys = ctx.nphys;
nf    = ctx.nf;
ny    = ctx.ny;
nx    = ctx.nx;
bsize = ctx.bsize;
phys  = ctx.phys_modes;
N     = nphys * bsize;

% Real -> complex physical state
q_complex = q_real(1:N) + 1i * q_real(N+1:2*N);

% Expand to full (nf x ny x nx) fields, filling in the conjugate modes
nvar = 4; if isfield(ctx,'nvar') && ~isempty(ctx.nvar); nvar = ctx.nvar; end
tmp.u = zeros(nf, ny, nx);
tmp.v = zeros(nf, ny, nx);
tmp.w = zeros(nf, ny, nx);
tmp.betavec  = ctx.betavec;
tmp.omegavec = ctx.omegavec;
if nvar == 5
    tmp.rho = zeros(nf, ny, nx);
    tmp.T   = zeros(nf, ny, nx);
    tmp.p   = zeros(nf, ny, nx);
    rho_b = ctx.StabGrid.rho;  T_b = ctx.StabGrid.T;   % base (ny x nx)
    gM2   = ctx.StabGrid.gamma * ctx.StabGrid.Ma^2;    % linear p' = (rho_b*T'+rho'*T_b)/gM2
end

for k = 1:nphys
    j   = phys(k);
    jc  = nf + 1 - j;
    off = (k-1)*bsize;
    qk_mat = reshape(q_complex(off+1:off+bsize), nvar*ny, nx);
    uj = qk_mat(1:ny,        :);
    vj = qk_mat(ny+1:2*ny,   :);
    wj = qk_mat(2*ny+1:3*ny, :);
    tmp.u(j,:,:) = uj;  tmp.v(j,:,:) = vj;  tmp.w(j,:,:) = wj;
    if nvar == 5
        rhoj = qk_mat(3*ny+1:4*ny, :);
        Tj   = qk_mat(4*ny+1:5*ny, :);
        pj   = (rho_b.*Tj + rhoj.*T_b) / gM2;          % linear pressure reconstruction
        tmp.rho(j,:,:) = rhoj;  tmp.T(j,:,:) = Tj;  tmp.p(j,:,:) = pj;
    end
    if jc ~= j                      % jc == j only for the (0,0) mean-flow mode
        tmp.u(jc,:,:) = conj(uj);
        tmp.v(jc,:,:) = conj(vj);
        tmp.w(jc,:,:) = conj(wj);
        if nvar == 5
            tmp.rho(jc,:,:) = conj(rhoj);
            tmp.T(jc,:,:)   = conj(Tj);
            tmp.p(jc,:,:)   = conj(pj);
        end
    end
end

%% NLT forcing (computed only for physical modes inside NLT_HNS)
f_all = NLT_HNS(ctx.StabGrid, 1:nf, tmp, ctx.D1, ctx.HB);

% Damp the NLT forcing through the outflow buffer, matching what the Picard path
% does, so the two solvers see the same forcing near the outflow.
ib = ctx.Bound.ibuff;
if ib <= nx
    bufcf_vec = repelem(ctx.Bound.bufcf(ib:nx), nvar*ny);
    idx = (ib-1)*nvar*ny+1 : nx*nvar*ny;
    f_all(:, idx) = f_all(:, idx) .* bufcf_vec;
end

%% Assemble residual
% Per-physical-mode residual (complex), with continuation factor AF on BC.
AF = ctx.AF;
F_complex = zeros(N, 1);
for k = 1:nphys
    j   = phys(k);
    off = (k-1)*bsize;
    qk  = q_complex(off+1:off+bsize);
    fk  = f_all(j,:).';
    rBC = AF * ctx.rBC{k};
    fk(rBC ~= 0) = 0;               % BC rows are Dirichlet: no NLT forcing there
    % F = ML*q - BC - NLT, the paper's M*Q - b(A0) - r. lambda enters only through
    % rBC, which is why the residual is affine in it and F_lambda is differenced once.
    F_complex(off+1:off+bsize) = ctx.ML{k}*qk - rBC - fk;
end

% Complex -> real
F_real = [real(F_complex); imag(F_complex)];
end
