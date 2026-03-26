function [MLAMFD,MLA1,MLA2] = LHS(BF,Bound,Grid,Init,Stab,StabRes,Opt)
%LHS  Assemble the global sparse LHS matrices for the HNS system.
%
%   [MLAMFD, MLA1, MLA2] = LHS(BF, Bound, Grid, Init, Stab, StabRes)
%
%   Constructs the block-pentadiagonal LHS operator from local coefficient
%   matrices (Ai..Gi) at each streamwise station using 4th-order FD stencils.
%   Sparse triplets are accumulated in cell arrays (O(n) cost).
%
%   Outputs:
%     MLAMFD  [N x N] sparse — LHS for the mean-flow-distortion (omega=0, beta=0)
%     MLA1    [N x N] sparse — base LHS (common to all modes)
%     MLA2    [N x N] sparse — omega-dependent part (ML = MLA1 + omega*MLA2)
%   where N = 4 * ny * nx.

% Pre-allocate cell arrays for sparse triplets (avoids O(n^2) concatenation)
nx = Grid.nx;
ny = Grid.ny;
bsize = 4*ny;

cell_i1MFD = cell(nx,1);  cell_j1MFD = cell(nx,1);  cell_s1MFD = cell(nx,1);
cell_i1    = cell(nx,1);  cell_j1    = cell(nx,1);  cell_s1    = cell(nx,1);
cell_i2    = cell(nx,1);  cell_j2    = cell(nx,1);  cell_s2    = cell(nx,1);
cell_i3    = cell(nx,1);  cell_j3    = cell(nx,1);  cell_s3    = cell(nx,1);
cell_i4    = cell(nx,1);  cell_j4    = cell(nx,1);  cell_s4    = cell(nx,1);

bufc   = Bound.bufc;
bufp   = Bound.bufp;
bufs   = Bound.bufs;
bufsp  = Bound.bufsp;

iu = Init.iu;
D1 = Grid.D1;
D2 = Grid.D2;
I  = Stab.I;
Z  = Stab.Z;
dxi = Grid.dxi;
Re = BF.Re;

% Pre-allocate local block matrices
% Chebyshev D1/D2 are dense → full arithmetic is faster for large ny
% FD D1/D2 are sparse (bandwidth 5) → keep sparse to exploit structure
D1_is_sparse = issparse(Grid.D1) && (nnz(Grid.D1) < 0.2 * ny^2);
use_full = (ny >= 60) && ~D1_is_sparse;

ML1MFD = spalloc(bsize, bsize*nx, 4*(3*ny^2+4*ny)+7*ny^2+4*ny);
ML1    = spalloc(bsize, bsize*nx, 4*(3*ny^2+4*ny)+7*ny^2+4*ny);
ML2    = spalloc(bsize, bsize*nx, 3*ny);
ML3    = spalloc(bsize, bsize*nx, 5*ny);
ML4    = spalloc(bsize, bsize*nx, 3*ny);

% Convert operators to full if beneficial
if use_full
    D1f = full(D1);  D2f = full(D2);  If = full(I);  Zf = full(Z);
else
    D1f = D1;  D2f = D2;  If = I;  Zf = Z;
end

% Precompute BC row indices.
% Chebyshev/Malik grid: index 1 = freestream (y=H), index ny = wall (y=0).
%   freestream row for variable k:  (k-1)*ny+1  (k=1,2,3 → rows 1, ny+1, 2*ny+1)
%   wall row for variable k:         k*ny       (k=1,2,3 → rows ny, 2*ny, 3*ny)
% For MFD (omega=0,beta=0), v-freestream BC is transpiration → skip k=2 freestream row.
bc_top     = [1, ny+1, 2*ny+1];   % freestream rows: all three velocity components
bc_top_mfd = [1, 2*ny+1];         % freestream rows: MFD (skip v, k=2)
bc_wall    = [ny, 2*ny, 3*ny];    % wall rows: all three velocity components
neumann_top = isfield(Opt,'bc_top') && strcmpi(Opt.bc_top,'neumann');

% --- Before the loop ---
t_lhs_start = tic;
lastPerc = 0;

for i=1:nx

 %% Prepare LHS building — use diag() for full, spdiags() for sparse

    if use_full
        U   = diag(BF.U(:,i));
        dxU = diag(BF.dxU(:,i));
        dyU = diag(BF.dyU(:,i));
        V   = diag(BF.V(:,i));
        dxV = diag(BF.dxV(:,i));
        dyV = diag(BF.dyV(:,i));
        W   = diag(BF.W(:,i));
        dxW = diag(BF.dxW(:,i));
        dyW = diag(BF.dyW(:,i));
    else
        U   = spdiags(BF.U(:,i),   0, ny, ny);
        dxU = spdiags(BF.dxU(:,i), 0, ny, ny);
        dyU = spdiags(BF.dyU(:,i), 0, ny, ny);
        V   = spdiags(BF.V(:,i),   0, ny, ny);
        dxV = spdiags(BF.dxV(:,i), 0, ny, ny);
        dyV = spdiags(BF.dyV(:,i), 0, ny, ny);
        W   = spdiags(BF.W(:,i),   0, ny, ny);
        dxW = spdiags(BF.dxW(:,i), 0, ny, ny);
        dyW = spdiags(BF.dyW(:,i), 0, ny, ny);
    end

 % Local transformations
    xix  = Grid.xix(:,i);
    etax = Grid.etax(:,i);
    xiy  = Grid.xiy(:,i);
    etay = Grid.etay(:,i);
    xixx = Grid.xixx(:,i);
    xiyy = Grid.xiyy(:,i);
    etaxx= Grid.etaxx(:,i);
    etayy= Grid.etayy(:,i);

 % Buffer
    bufs_i  = bufs(:,i);
    bufsp_i = bufsp(:,i);
    bufc_i  = bufc(i);
    bufp_i  = bufp(i);

 % Create common terms (ct) for matrices A-F
    Act =  U.*etax*D1f...
          +V.*etay*D1f...
          -1/Re.*(etax.^2+etay.^2).*D2f...
          -1/Re.*(etaxx+etayy).*D1f;

    Ect =  U.*xix...
          +V.*xiy...
          -1/Re*(xixx+xiyy).*I...
          -1/Re*(2*etax.*xix+2*etay.*xiy).*D1f;

%% construct coefficient matrices
    Ai=[(Act+dxU).*bufs_i       dyU.*bufs_i        Zf    etax.*D1f.*bufs_i
            dxV.*bufs_i        (Act+dyV).*bufs_i   Zf    etay.*D1f.*bufs_i
            dxW.*bufs_i         dyW.*bufs_i       Act.*bufs_i Zf
           (etax.*D1f).*bufsp_i (etay.*D1f).*bufsp_i Zf             Zf];

    Bi=[-iu*If.*bufs_i  Zf  Zf  Zf
         Zf  -iu*If.*bufs_i Zf  Zf
         Zf  Zf  -iu*If.*bufs_i Zf
         Zf  Zf                Zf  Zf];

    Ci=[iu*W.*bufs_i    Zf   Zf    Zf
         Zf    iu*W.*bufs_i  Zf    Zf
         Zf    Zf   iu*W.*bufs_i   iu*If.*bufs_i
         Zf    Zf   iu*If.*bufsp_i  Zf];

    Di=[If/Re.*bufs_i  Zf  Zf  Zf
         Zf  If/Re.*bufs_i Zf  Zf
         Zf  Zf  If/Re.*bufs_i Zf
         Zf                      Zf  Zf  Zf];

    Ei=[Ect.*bufs_i      Zf        Zf    If.*xix.*bufs_i
        Zf       Ect.*bufs_i       Zf    If.*xiy.*bufs_i
        Zf       Zf        Ect.*bufs_i            Zf
        If.*xix.*bufsp_i   If.*xiy.*bufsp_i Zf  Zf];

    Fi=[-If/Re.*(xix.^2+xiy.^2).*bufs_i.*bufp_i   Zf       Zf       Zf
        Zf       -If/Re.*(xix.^2+xiy.^2).*bufs_i.*bufp_i   Zf       Zf
        Zf       Zf       -If/Re.*(xix.^2+xiy.^2).*bufs_i.*bufp_i   Zf
        Zf       Zf                                                   Zf       Zf];

    % Ensure pressure is copied downward from the surface
    if use_full
        Pres = diag(1-bufsp_i) + diag(-1*(1-bufsp(2:end,i)), -1);
    else
        Pres = I.*(1-bufsp_i) + spdiags(-1*(1-bufsp(2:end,i)), -1, ny, ny);
    end

    % Set up buffer and embedded boundary method
    Gi=[If.*(1-bufs_i)+If.*bufs(:,end)*(1-bufc_i) Zf Zf Zf
        Zf If.*(1-bufs_i) Zf Zf
        Zf Zf If.*(1-bufs_i)+If.*bufs(:,end)*(1-bufc_i)  Zf
        Zf Zf Zf Pres];



    % construct the ML block for station i

    sten=bsize*(i-1)+1:bsize*i; % row coordinates of FD stencil
    sten1=1:bsize;              % stencil height
    bc_col0 = bsize*(i-1);     % column offset for BC diagonal entries

    if i==1
        ML1MFD(sten1,sten)=eye(bsize);
        ML1(sten1,sten)=eye(bsize);
        ML2(sten1,sten)=zeros(bsize);
        ML3(sten1,sten)=zeros(bsize);
        ML4(sten1,sten)=zeros(bsize);

    elseif (i==2)

        stenr=[sten-bsize sten sten+bsize sten+2*bsize sten+3*bsize];

        ML1MFD(sten1,stenr)=[-3*Ei/(12*dxi)+11*Fi/(12*dxi^2), ...
                         (Ai-10*Ei/(12*dxi)-20*Fi/(12*dxi^2))+Gi/bufc_i,...
                             18*Ei/(12*dxi)+ 6*Fi/(12*dxi^2), ...
                            - 6*Ei/(12*dxi)+ 4*Fi/(12*dxi^2), ...
                                Ei/(12*dxi)- 1*Fi/(12*dxi^2)]*bufc_i;

        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi*bufc_i;
        ML3(sten1,sten)=Ci*bufc_i;
        ML4(sten1,sten)=Di*bufc_i;

        ML1MFD(bc_top_mfd,:)=0; ML1(bc_top,:)=0; ML2(bc_top,:)=0; ML3(bc_top,:)=0; ML4(bc_top,:)=0;
        if neumann_top
            for k=1:3
                r = (k-1)*ny+1;
                cs = bc_col0+(k-1)*ny+1;  ce = bc_col0+k*ny;
                ML1(r,cs:ce) = D1f(1,:);
                if k~=2; ML1MFD(r,cs:ce) = D1f(1,:); end  % skip v-freestream for MFD (transpiration)
            end
        else
            for r=bc_top_mfd; ML1MFD(r,bc_col0+r)=1; end
            for r=bc_top;     ML1(r,   bc_col0+r)=1; end
        end
        ML1MFD(bc_wall,:)=0; ML1(bc_wall,:)=0; ML2(bc_wall,:)=0; ML3(bc_wall,:)=0; ML4(bc_wall,:)=0;
        for r=bc_wall; ML1MFD(r,bc_col0+r)=1; ML1(r,bc_col0+r)=1; end

    elseif (i==nx-1)

        stenr=[sten-2*bsize sten-bsize sten sten+bsize];

        ML1MFD(sten1,stenr)=[Ei/(6*dxi), ...
                          -6*Ei/(6*dxi)+  Fi/(dxi^2),...
                       (Ai+3*Ei/(6*dxi)-2*Fi/(dxi^2))+Gi/bufc_i,...
                           2*Ei/(6*dxi)+  Fi/(dxi^2)]*bufc_i;

        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi*bufc_i;
        ML3(sten1,sten)=Ci*bufc_i;
        ML4(sten1,sten)=Di*bufc_i;

        ML1MFD(bc_top_mfd,:)=0; ML1(bc_top,:)=0; ML2(bc_top,:)=0; ML3(bc_top,:)=0; ML4(bc_top,:)=0;
        if neumann_top
            for k=1:3
                r = (k-1)*ny+1;
                cs = bc_col0+(k-1)*ny+1;  ce = bc_col0+k*ny;
                ML1(r,cs:ce) = D1f(1,:);
                if k~=2; ML1MFD(r,cs:ce) = D1f(1,:); end  % skip v-freestream for MFD (transpiration)
            end
        else
            for r=bc_top_mfd; ML1MFD(r,bc_col0+r)=1; end
            for r=bc_top;     ML1(r,   bc_col0+r)=1; end
        end
        ML1MFD(bc_wall,:)=0; ML1(bc_wall,:)=0; ML2(bc_wall,:)=0; ML3(bc_wall,:)=0; ML4(bc_wall,:)=0;
        for r=bc_wall; ML1MFD(r,bc_col0+r)=1; ML1(r,bc_col0+r)=1; end

    elseif (i==nx)
        stenr=[sten-2*bsize sten-bsize sten];

        ML1MFD(sten1,stenr)=[Ei/(2*dxi)+  Fi/(dxi^2), ...
                          -4*Ei/(2*dxi)-2*Fi/(dxi^2), ...
                       (Ai+3*Ei/(2*dxi)+  Fi/(dxi^2))+Gi/bufc_i]*bufc_i;

        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi*bufc_i;
        ML3(sten1,sten)=Ci*bufc_i;
        ML4(sten1,sten)=Di*bufc_i;

        ML1MFD(bc_top_mfd,:)=0; ML1(bc_top,:)=0; ML2(bc_top,:)=0; ML3(bc_top,:)=0; ML4(bc_top,:)=0;
        if neumann_top
            for k=1:3
                r = (k-1)*ny+1;
                cs = bc_col0+(k-1)*ny+1;  ce = bc_col0+k*ny;
                ML1(r,cs:ce) = D1f(1,:);
                if k~=2; ML1MFD(r,cs:ce) = D1f(1,:); end  % skip v-freestream for MFD (transpiration)
            end
        else
            for r=bc_top_mfd; ML1MFD(r,bc_col0+r)=1; end
            for r=bc_top;     ML1(r,   bc_col0+r)=1; end
        end
        ML1MFD(bc_wall,:)=0; ML1(bc_wall,:)=0; ML2(bc_wall,:)=0; ML3(bc_wall,:)=0; ML4(bc_wall,:)=0;
        for r=bc_wall; ML1MFD(r,bc_col0+r)=1; ML1(r,bc_col0+r)=1; end
    elseif i == Bound.stepi-1 && Grid.StepH ~= 0
        yi = Bound.yi;

        % Values under step height
        sten1u=[0*ny+((yi+1):ny) 1*ny+((yi+1):ny) 2*ny+((yi+1):ny) 3*ny+((yi+1):ny) ];
        stenr=[sten-3*bsize sten-2*bsize sten-bsize sten sten+bsize];

        ML1MFD(sten1u,stenr)=[- 1*Ei(sten1u,:)/(12*dxi)-   Fi(sten1u,:)/(12*dxi^2),...
                               6*Ei(sten1u,:)/(12*dxi)+ 4*Fi(sten1u,:)/(12*dxi^2),...
                             -18*Ei(sten1u,:)/(12*dxi)+ 6*Fi(sten1u,:)/(12*dxi^2),...
                  Ai(sten1u,:)+10*Ei(sten1u,:)/(12*dxi)-20*Fi(sten1u,:)/(12*dxi^2)+Gi(sten1u,:)/bufc_i,...
                               3*Ei(sten1u,:)/(12*dxi)+11*Fi(sten1u,:)/(12*dxi^2)]*bufc_i;

        ML1(sten1u,stenr)=ML1MFD(sten1u,stenr);
        ML2(sten1u,sten)=Bi(sten1u,:)*bufc_i;
        ML3(sten1u,sten)=Ci(sten1u,:)*bufc_i;
        ML4(sten1u,sten)=Di(sten1u,:)*bufc_i;

        % Values above step
        sten1a=[0*ny+(1:yi) 1*ny+(1:yi) 2*ny+(1:yi) 3*ny+(1:yi) ];
        stenr=[sten-2*bsize sten-bsize sten sten+bsize sten+2*bsize];

        ML1MFD(sten1a,stenr)=[   Ei(sten1a,:)/(12*dxi)-   Fi(sten1a,:)/(12*dxi^2),...
                             -8*Ei(sten1a,:)/(12*dxi)+16*Fi(sten1a,:)/(12*dxi^2),...
                (Ai(sten1a,:)                        -30*Fi(sten1a,:)/(12*dxi^2))+Gi(sten1a,:)/bufc_i,...
                              8*Ei(sten1a,:)/(12*dxi)+16*Fi(sten1a,:)/(12*dxi^2),...
                             -  Ei(sten1a,:)/(12*dxi)-   Fi(sten1a,:)/(12*dxi^2)]*bufc_i;

        ML1(sten1a,stenr)=ML1MFD(sten1a,stenr);
        ML2(sten1a,sten)=Bi(sten1a,:)*bufc_i;
        ML3(sten1a,sten)=Ci(sten1a,:)*bufc_i;
        ML4(sten1a,sten)=Di(sten1a,:)*bufc_i;

        ML1MFD(bc_top_mfd,:)=0; ML1(bc_top,:)=0; ML2(bc_top,:)=0; ML3(bc_top,:)=0; ML4(bc_top,:)=0;
        if neumann_top
            for k=1:3
                r = (k-1)*ny+1;
                cs = bc_col0+(k-1)*ny+1;  ce = bc_col0+k*ny;
                ML1(r,cs:ce) = D1f(1,:);
                if k~=2; ML1MFD(r,cs:ce) = D1f(1,:); end  % skip v-freestream for MFD (transpiration)
            end
        else
            for r=bc_top_mfd; ML1MFD(r,bc_col0+r)=1; end
            for r=bc_top;     ML1(r,   bc_col0+r)=1; end
        end
        ML1MFD(bc_wall,:)=0; ML1(bc_wall,:)=0; ML2(bc_wall,:)=0; ML3(bc_wall,:)=0; ML4(bc_wall,:)=0;
        for r=bc_wall; ML1MFD(r,bc_col0+r)=1; ML1(r,bc_col0+r)=1; end

    elseif i == Bound.stepi && Grid.StepH ~= 0
        yi = Bound.yi;

        % Values under step height
        sten1u=[0*ny+((yi+1):ny) 1*ny+((yi+1):ny) 2*ny+((yi+1):ny) 3*ny+((yi+1):ny) ];
        stenr=[sten-4*bsize sten-3*bsize sten-2*bsize sten-bsize sten];

        ML1MFD(sten1u,stenr)=[+ 3*Ei(sten1u,:)/(12*dxi)+ 11*Fi(sten1u,:)/(12*dxi^2),...
                             -16*Ei(sten1u,:)/(12*dxi)- 56*Fi(sten1u,:)/(12*dxi^2),...
                             +36*Ei(sten1u,:)/(12*dxi)+114*Fi(sten1u,:)/(12*dxi^2),...
                             -48*Ei(sten1u,:)/(12*dxi)-104*Fi(sten1u,:)/(12*dxi^2),...
                 (Ai(sten1u,:)+25*Ei(sten1u,:)/(12*dxi)+ 35*Fi(sten1u,:)/(12*dxi^2))+Gi(sten1u,:)/bufc_i]*bufc_i;

        ML1(sten1u,stenr)=ML1MFD(sten1u,stenr);
        ML2(sten1u,sten)=Bi(sten1u,:)*bufc_i;
        ML3(sten1u,sten)=Ci(sten1u,:)*bufc_i;
        ML4(sten1u,sten)=Di(sten1u,:)*bufc_i;

        % Values above step
        sten1a=[0*ny+(1:yi) 1*ny+(1:yi) 2*ny+(1:yi) 3*ny+(1:yi) ];
        stenr=[sten-2*bsize sten-bsize sten sten+bsize sten+2*bsize];

        ML1MFD(sten1a,stenr)=[Ei(sten1a,:)/(12*dxi)- 1*Fi(sten1a,:)/(12*dxi^2), ...
                          -8*Ei(sten1a,:)/(12*dxi)+16*Fi(sten1a,:)/(12*dxi^2), ...
             (Ai(sten1a,:)                        -30*Fi(sten1a,:)/(12*dxi^2))+Gi(sten1a,:)/bufc_i, ...
                           8*Ei(sten1a,:)/(12*dxi)+16*Fi(sten1a,:)/(12*dxi^2), ...
                          -1*Ei(sten1a,:)/(12*dxi)- 1*Fi(sten1a,:)/(12*dxi^2)]*bufc_i;

        ML1(sten1a,stenr)=ML1MFD(sten1a,stenr);
        ML2(sten1a,sten)=Bi(sten1a,:)*bufc_i;
        ML3(sten1a,sten)=Ci(sten1a,:)*bufc_i;
        ML4(sten1a,sten)=Di(sten1a,:)*bufc_i;

        ML1MFD(bc_top_mfd,:)=0; ML1(bc_top,:)=0; ML2(bc_top,:)=0; ML3(bc_top,:)=0; ML4(bc_top,:)=0;
        if neumann_top
            for k=1:3
                r = (k-1)*ny+1;
                cs = bc_col0+(k-1)*ny+1;  ce = bc_col0+k*ny;
                ML1(r,cs:ce) = D1f(1,:);
                if k~=2; ML1MFD(r,cs:ce) = D1f(1,:); end  % skip v-freestream for MFD (transpiration)
            end
        else
            for r=bc_top_mfd; ML1MFD(r,bc_col0+r)=1; end
            for r=bc_top;     ML1(r,   bc_col0+r)=1; end
        end
        ML1MFD(bc_wall,:)=0; ML1(bc_wall,:)=0; ML2(bc_wall,:)=0; ML3(bc_wall,:)=0; ML4(bc_wall,:)=0;
        for r=bc_wall; ML1MFD(r,bc_col0+r)=1; ML1(r,bc_col0+r)=1; end
    else
        stenr=[sten-2*bsize sten-bsize sten sten+bsize sten+2*bsize];

        ML1MFD(sten1,stenr)=[Ei/(12*dxi)- 1*Fi/(12*dxi^2), ...
                          -8*Ei/(12*dxi)+16*Fi/(12*dxi^2), ...
                        (Ai             -30*Fi/(12*dxi^2))+Gi/bufc_i, ...
                           8*Ei/(12*dxi)+16*Fi/(12*dxi^2), ...
                            -Ei/(12*dxi)- 1*Fi/(12*dxi^2)]*bufc_i;

        ML1(sten1,stenr)=ML1MFD(sten1,stenr);
        ML2(sten1,sten)=Bi*bufc_i;
        ML3(sten1,sten)=Ci*bufc_i;
        ML4(sten1,sten)=Di*bufc_i;

        ML1MFD(bc_top_mfd,:)=0; ML1(bc_top,:)=0; ML2(bc_top,:)=0; ML3(bc_top,:)=0; ML4(bc_top,:)=0;
        if neumann_top
            for k=1:3
                r = (k-1)*ny+1;
                cs = bc_col0+(k-1)*ny+1;  ce = bc_col0+k*ny;
                ML1(r,cs:ce) = D1f(1,:);
                if k~=2; ML1MFD(r,cs:ce) = D1f(1,:); end  % skip v-freestream for MFD (transpiration)
            end
        else
            for r=bc_top_mfd; ML1MFD(r,bc_col0+r)=1; end
            for r=bc_top;     ML1(r,   bc_col0+r)=1; end
        end
        ML1MFD(bc_wall,:)=0; ML1(bc_wall,:)=0; ML2(bc_wall,:)=0; ML3(bc_wall,:)=0; ML4(bc_wall,:)=0;
        for r=bc_wall; ML1MFD(r,bc_col0+r)=1; ML1(r,bc_col0+r)=1; end
    end

% --- Progress ---
perc = floor(i * 100 / nx);
if perc >= lastPerc + 25 && perc < 100
    elapsed = toc(t_lhs_start);
    eta = elapsed / i * (nx - i);
    fprintf('   LHS: %d%% (%.0fs elapsed, ~%.0fs remaining)\n', perc, elapsed, eta);
    lastPerc = perc;
end

% Store triplets in cell arrays (O(1) per iteration instead of O(n))
[ii,jj,ss] = find(ML1MFD); cell_i1MFD{i} = ii+bsize*(i-1); cell_j1MFD{i} = jj; cell_s1MFD{i} = ss;
[ii,jj,ss] = find(ML1);    cell_i1{i}    = ii+bsize*(i-1); cell_j1{i}    = jj; cell_s1{i}    = ss;
[ii,jj,ss] = find(ML2);    cell_i2{i}    = ii+bsize*(i-1); cell_j2{i}    = jj; cell_s2{i}    = ss;
[ii,jj,ss] = find(ML3);    cell_i3{i}    = ii+bsize*(i-1); cell_j3{i}    = jj; cell_s3{i}    = ss;
[ii,jj,ss] = find(ML4);    cell_i4{i}    = ii+bsize*(i-1); cell_j4{i}    = jj; cell_s4{i}    = ss;

% Reset local blocks
ML1MFD(:) = 0; ML1(:) = 0; ML2(:) = 0; ML3(:) = 0; ML4(:) = 0;

end

% --- After the loop ---
fprintf('   LHS: 100%% (%.1fs total)\n', toc(t_lhs_start));

% Concatenate cell arrays and clear immediately to reduce peak memory
N = bsize * nx;

MLAMFD = sparse(cell2mat(cell_i1MFD), cell2mat(cell_j1MFD), cell2mat(cell_s1MFD), N, N);
clear cell_i1MFD cell_j1MFD cell_s1MFD

MLA1 = sparse(cell2mat(cell_i1), cell2mat(cell_j1), cell2mat(cell_s1), N, N);
clear cell_i1 cell_j1 cell_s1

if any(StabRes.omegavec)
    MLA2 = sparse(cell2mat(cell_i2), cell2mat(cell_j2), cell2mat(cell_s2), N, N);
end
clear cell_i2 cell_j2 cell_s2

if any(StabRes.betavec)
    MLA3 = sparse(cell2mat(cell_i3), cell2mat(cell_j3), cell2mat(cell_s3), N, N);
    MLA4 = sparse(cell2mat(cell_i4), cell2mat(cell_j4), cell2mat(cell_s4), N, N);
end
clear cell_i3 cell_j3 cell_s3 cell_i4 cell_j4 cell_s4
clear ML1MFD ML1 ML2 ML3 ML4

end
