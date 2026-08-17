function varargout = block_thomas(op, varargin)
%BLOCK_THOMAS  Dual-Schur block-pentadiagonal Thomas solver (Paredes et al. 2024).
%   Direct solver for the station-major block-pentadiagonal HNS operator ML
%   (block size b = Nq*Ny; a block is nonzero only for |i-j|<=2, a consequence
%   of the 4th-order streamwise stencil). The domain is split into np
%   partitions separated by two-station interfaces. Each partition is
%   factored with an unpivoted block-Thomas recursion (the dense diagonal
%   blocks are themselves inverted with pivoted LAPACK), and the partitions
%   are coupled through a Schur complement on the interfaces. The
%   across-station recursion is unpivoted but stays stable because each
%   partition only marches N/np stations, which bounds the parasitic growth
%   that an elliptic operator picks up under a forward sweep (a single sweep,
%   np=1, is NOT stable). Default np=2 (fewest partitions); verified exact
%   against the banded LU factorization through nx=1600.
%
%   Operations (the first argument selects the operation):
%     x  = block_thomas('solve',  ML, b, rhs, np, oocdir)   one-shot factor+apply+free
%     F  = block_thomas('factor', ML, b, np, oocdir)        build a reusable factor
%     x  = block_thomas('apply',  F, rhs)                   solve ML*x = rhs from F
%          block_thomas('free',   F)                        delete an F's on-disk files
%          block_thomas('free',   OOCDIR, 'reg')            register an OOCDIR for a sweep
%     mk = block_thomas('free',   [], 'mark')               registry checkpoint
%          block_thomas('free',   mk, 'sweep')              delete dirs registered after mk
%
%   Pass OOCDIR to stream the factor to disk (out-of-core, OOC): only a few
%   blocks are ever held in RAM at once (Paredes et al. 2024), so resolutions
%   whose in-core factor would exceed memory stay reachable. rhs may have
%   multiple columns. The registry lets main() install a cleanup routine
%   (onCleanup) that removes a run's streamed files whenever MATLAB exits.

switch lower(op)
    case 'solve';  [varargout{1:nargout}] = bt_solve(varargin{:});
    case 'factor'; [varargout{1:nargout}] = bt_factor(varargin{:});
    case 'apply';  [varargout{1:nargout}] = bt_apply(varargin{:});
    case 'free';   [varargout{1:nargout}] = bt_free(varargin{:});
    otherwise
        error('block_thomas:unknownOp', ...
            'unknown operation ''%s'' (use solve/factor/apply/free)', op);
end
end

%% One-shot solve
function x = bt_solve(ML, b, rhs, np, oocdir)
if nargin < 4 || isempty(np); np = 2; end
if nargin < 5; oocdir = ''; end
F = bt_factor(ML, b, np, oocdir);
x = bt_apply(F, rhs);
bt_free(F);
end

%% Factor: build a right-hand-side (RHS) independent dual-Schur factor of ML (in-core, or streamed to OOCDIR)
function F = bt_factor(ML, b, np, oocdir)
if nargin < 3 || isempty(np); np = 2; end
if nargin < 4; oocdir = ''; end
ooc = ~isempty(oocdir);
if ooc && ~exist(oocdir,'dir'); mkdir(oocdir); end
n = size(ML,1); N = round(n/b);

% stations whose stencil reaches past the pentadiagonal band (near a boundary) go into the interface
[r,c] = find(ML); bi = ceil(r/b); bj = ceil(c/b); far = abs(bi-bj) > 2;
bad = unique([bi(far); bj(far)]);
istart = 1; iend = N;
if ~isempty(bad)
    lo = bad(bad <= N/2); if ~isempty(lo); istart = max(lo) + 1; end
    hi = bad(bad >  N/2); if ~isempty(hi); iend  = min(hi) - 1; end
end

% partition the interior stations; the boundary blocks plus the seams between partitions form the interface
label = zeros(N,1);
Nin = iend - istart + 1; nint = Nin - 2*(np-1); chunk = floor(nint/np); pos = istart;
for p = 1:np
    csz = chunk + double(p <= mod(nint,np));
    label(pos:pos+csz-1) = p; pos = pos + csz;
    if p < np; label(pos:pos+1) = 0; pos = pos + 2; end
end
sd = @(st) reshape(((st(:)-1)*b + (1:b)).', [], 1);
ifDof = sd(find(label==0));

% factor each partition, accumulate the interface Schur complement
St = ML(ifDof,ifDof);
Cp = cell(np,1); W = cell(np,1); wcol = cell(np,1); pdof = cell(np,1);
part = cell(np,1); Dp = cell(np,1); wfile = cell(np,1); wsize = zeros(np,2);
for p = 1:np
    pd = sd(find(label==p)); pdof{p} = pd;
    Dp{p} = ML(pd,pd);
    Bp = ML(pd,ifDof); Cp{p} = ML(ifDof,pd);
    cols = find(any(Bp,1)); wcol{p} = cols; ncols = numel(cols);
    Np = round(numel(pd)/b);
    if ooc
        facfile = fullfile(oocdir, sprintf('fac_p%03d.bin', p));
        wfile{p} = fullfile(oocdir, sprintf('w_p%03d.bin', p));
        part{p} = local_thomas_factor(Dp{p}, b, full(Bp(:,cols)), facfile, wfile{p});
        wsize(p,:) = [numel(pd), ncols];
        % accumulate St -= Cp*W by streaming W back from disk one block at a time
        fidW = fopen(wfile{p},'r'); recW = 2*b*ncols;
        for i = 1:Np
            fseek(fidW,(Np-i)*recW*8,'bof'); v = bt_fread(fidW,recW,wfile{p});
            qi = reshape(complex(v(1:b*ncols), v(b*ncols+1:end)), b, ncols);
            rows = (i-1)*b+(1:b);
            St(:,cols) = St(:,cols) - Cp{p}(:,rows)*qi;
        end
        fclose(fidW);
    else
        [W{p}, part{p}] = local_thomas_factor(Dp{p}, b, full(Bp(:,cols)));
        St(:,cols) = St(:,cols) - Cp{p}*W{p};
    end
end
[StL, StU, Stp] = lu(full(St), 'vector');

F = struct('b',b, 'n',n, 'np',np, 'ifDof',ifDof, 'pdof',{pdof}, 'wcol',{wcol}, ...
           'Cp',{Cp}, 'part',{part}, 'Dp',{Dp}, ...
           'StL',StL, 'StU',StU, 'Stp',Stp, 'ooc',ooc, 'oocdir',oocdir);
if ooc
    F.wfile = wfile; F.wsize = wsize;
else
    F.W = W;
end
end

function varargout = local_thomas_factor(D, b, Bcols, facfile, wfile)
%LOCAL_THOMAS_FACTOR  Factor one partition and produce its interface response W.
%   In-core (3 args): returns [W, part], part = {L_i,U_i,pv_i,Y_i,Z_i}.
%   Out-of-core (5 args): streams the station factor blocks [L,U,Y,Z] to FACFILE
%   and W block-by-block to WFILE (pv kept in RAM), holding only a few blocks in
%   memory at a time; returns part = {pv, N, facfile}. Nothing the size of a
%   whole partition is ever materialized in memory.
stream = nargin >= 4 && ~isempty(facfile);
n = size(D,1); N = round(n/b); zb = sparse(b,b); nc = size(Bcols,2);
if stream
    fid = fopen(facfile,'w'); Pc = cell(N,1);
    wfwd = [facfile '.wf']; fidWf = fopen(wfwd,'w');   % scratch file for the forward W sweep, deleted after use
else
    Lc = cell(N,1); Uc = cell(N,1); Pc = cell(N,1); Yc = cell(N,1); Zc = cell(N,1);
    Wc = cell(N,1);
end
Ym1=zb; Ym2=zb; Zm1=zb; Zm2=zb; rm1=zeros(b,nc); rm2=zeros(b,nc);
blk = @(i,j) D((i-1)*b+(1:b),(j-1)*b+(1:b));
for i = 1:N
    if i>=3, A = blk(i,i-2); else, A = zb; end
    if i>=2, Bb = blk(i,i-1); else, Bb = zb; end
    C = blk(i,i);
    if i<=N-1, Dd = blk(i,i+1); else, Dd = zb; end
    if i<=N-2, E  = blk(i,i+2); else, E  = zb; end
    K = Bb - A*Ym2;
    G = full(C - K*Ym1 - A*Zm2);
    [L,U,pv] = lu(G,'vector');
    solG = @(V) U\(L\V(pv,:));
    Yi = solG(full(Dd - K*Zm1));
    Zi = solG(full(E));
    Wi = solG(Bcols((i-1)*b+(1:b),:) - K*rm1 - A*rm2);
    if stream
        bt_fwrite(fid,[real(L(:));imag(L(:));real(U(:));imag(U(:)); ...
                    real(Yi(:));imag(Yi(:));real(Zi(:));imag(Zi(:))],facfile);
        bt_fwrite(fidWf,[real(Wi(:));imag(Wi(:))],wfwd);    % forward W to disk
        Pc{i} = pv;
    else
        Lc{i}=L; Uc{i}=U; Pc{i}=pv; Yc{i}=Yi; Zc{i}=Zi; Wc{i}=Wi;
    end
    Ym2=Ym1; Ym1=Yi; Zm2=Zm1; Zm1=Zi; rm2=rm1; rm1=Wi;
end
if stream; fclose(fid); fclose(fidWf); end

% back-substitution for the interface response W = D \ Bcols
qp1=zeros(b,nc); qp2=zeros(b,nc);
if stream
    fid   = fopen(facfile,'r'); rec = 8*b*b; bb = b*b;
    fidWf = fopen(wfwd,'r');    recWf = 2*b*nc;
    fidW  = fopen(wfile,'w');   recW = 2*b*nc;   % W is written out contiguously, in reverse station order
    for i = N:-1:1                               % record (N-i) in the file holds station i
        fseek(fid,(i-1)*rec*8,'bof'); v = bt_fread(fid,rec,facfile);
        Yi = reshape(complex(v(4*bb+1:5*bb), v(5*bb+1:6*bb)), b, b);
        Zi = reshape(complex(v(6*bb+1:7*bb), v(7*bb+1:8*bb)), b, b);
        fseek(fidWf,(i-1)*recWf*8,'bof'); vw = bt_fread(fidWf,recWf,wfwd);
        Wi = reshape(complex(vw(1:b*nc), vw(b*nc+1:end)), b, nc);
        qi = Wi - Yi*qp1 - Zi*qp2;
        bt_fwrite(fidW,[real(qi(:));imag(qi(:))],wfile);
        qp2=qp1; qp1=qi;
    end
    fclose(fid); fclose(fidWf); fclose(fidW); delete(wfwd);
    part = struct('pv',{Pc}, 'N',N, 'facfile',facfile);
    varargout = {part};
else
    W = zeros(n,nc);
    for i = N:-1:1
        qi = Wc{i} - Yc{i}*qp1 - Zc{i}*qp2;
        W((i-1)*b+(1:b),:) = qi; qp2=qp1; qp1=qi;
    end
    part = struct('L',{Lc}, 'U',{Uc}, 'pv',{Pc}, 'Y',{Yc}, 'Z',{Zc});
    varargout = {W, part};
end
end

%% Apply: solve ML*x = rhs from a precomputed factor F (rhs may have many columns)
function x = bt_apply(F, rhs)
b = F.b; n = F.n; np = F.np; nc = size(rhs,2);
ft = rhs(F.ifDof,:);
up = cell(np,1);
for p = 1:np
    up{p} = local_thomas_apply(F.part{p}, F.Dp{p}, b, rhs(F.pdof{p},:), F.ooc);
    ft = ft - F.Cp{p}*up{p};
end
qs = F.StU \ (F.StL \ ft(F.Stp,:));
x = zeros(n, nc); x(F.ifDof,:) = qs;
for p = 1:np
    if F.ooc
        ncols = F.wsize(p,2); Np = round(F.wsize(p,1)/b);
        qcols = qs(F.wcol{p},:); pd = F.pdof{p};
        fidW = fopen(F.wfile{p},'r'); recW = 2*b*ncols;
        for i = 1:Np
            fseek(fidW,(Np-i)*recW*8,'bof'); v = bt_fread(fidW,recW,F.wfile{p});
            qi = reshape(complex(v(1:b*ncols), v(b*ncols+1:end)), b, ncols);
            rows = (i-1)*b+(1:b);
            x(pd(rows),:) = up{p}(rows,:) - qi*qcols;
        end
        fclose(fidW);
    else
        x(F.pdof{p},:) = up{p} - F.W{p}*qs(F.wcol{p},:);
    end
end
end

function up = local_thomas_apply(part, D, b, rhsp, ooc)
%LOCAL_THOMAS_APPLY  Compute the particular solution up = D \ rhsp on one partition.
%   K_i and A_i are recomputed from the sparse partition matrix D and a rolling
%   Y_{i-2} buffer, so only the G_i LU factors and Y_i, Z_i need to be kept
%   around (read from disk when OOC).
n = size(rhsp,1); N = round(n/b); zb = sparse(b,b); nc = size(rhsp,2);
rm1 = zeros(b,nc); rm2 = zeros(b,nc); rr = cell(N,1);
blk = @(i,j) D((i-1)*b+(1:b),(j-1)*b+(1:b));
Yb1 = zb; Yb2 = zb;                      % Y_{i-1}, Y_{i-2}
if ooc; fid = fopen(part.facfile,'r'); rec = 8*b*b; bb = b*b; end
for i = 1:N
    if i>=3, A = blk(i,i-2); Ym2 = Yb2; else, A = zb; Ym2 = zb; end
    if i>=2, Bb = blk(i,i-1); else, Bb = zb; end
    K = Bb - A*Ym2;
    if ooc
        fseek(fid,(i-1)*rec*8,'bof'); v = bt_fread(fid,4*bb,part.facfile);   % L,U
        L  = reshape(complex(v(1:bb),      v(bb+1:2*bb)),   b, b);
        U  = reshape(complex(v(2*bb+1:3*bb), v(3*bb+1:4*bb)), b, b);
        fseek(fid,((i-1)*rec+4*bb)*8,'bof'); vy = bt_fread(fid,2*bb,part.facfile); % Y
        Yi = reshape(complex(vy(1:bb), vy(bb+1:2*bb)), b, b);
        pv = part.pv{i};
    else
        L = part.L{i}; U = part.U{i}; pv = part.pv{i}; Yi = part.Y{i};
    end
    V = rhsp((i-1)*b+(1:b),:) - K*rm1 - A*rm2;
    rr{i} = U \ (L \ V(pv,:));
    Yb2 = Yb1; Yb1 = Yi; rm2 = rm1; rm1 = rr{i};
end
up = zeros(n,nc); qp1 = zeros(b,nc); qp2 = zeros(b,nc);
for i = N:-1:1
    if ooc
        fseek(fid,((i-1)*rec+4*bb)*8,'bof'); v = bt_fread(fid,4*bb,part.facfile); % Y,Z
        Yi = reshape(complex(v(1:bb),      v(bb+1:2*bb)),   b, b);
        Zi = reshape(complex(v(2*bb+1:3*bb), v(3*bb+1:4*bb)), b, b);
    else
        Yi = part.Y{i}; Zi = part.Z{i};
    end
    qi = rr{i} - Yi*qp1 - Zi*qp2;
    up((i-1)*b+(1:b),:) = qi; qp2 = qp1; qp1 = qi;
end
if ooc; fclose(fid); end
end

%% Free / out-of-core temporary-file registry (see the header above for what each mode does)
function out = bt_free(F, mode)
persistent reg
if isempty(reg); reg = {}; end
out = [];
if nargin >= 2
    switch mode
        case 'reg'
            reg{end+1} = F;                          % F is an oocdir string
            return
        case 'mark'
            out = numel(reg);
            return
        case {'sweep','sweep_all'}
            mk = 0;
            if strcmp(mode,'sweep') && ~isempty(F); mk = F; end
            for i = numel(reg):-1:mk+1
                if ~isempty(reg{i}) && exist(reg{i},'dir'); rmdir(reg{i},'s'); end
                reg(i) = [];
            end
            return
    end
end

% default: free a factor struct's files
if isstruct(F) && isfield(F,'ooc') && F.ooc && isfield(F,'oocdir') ...
        && ~isempty(F.oocdir) && exist(F.oocdir,'dir')
    rmdir(F.oocdir,'s');
end
end

%% Guarded out-of-core I/O (turns a disk-full truncation into a clear error message)
function v = bt_fread(fid, n, fname)
v = fread(fid, n, 'double');
if numel(v) ~= n
    error('block_thomas:oocTruncated', ...
        ['out-of-core read from ''%s'' returned %d of %d values: the file is ' ...
         'truncated, most likely from a full disk during factorization. Free disk ' ...
         'space, lower the resolution, or use lu_mode=''lapack_band''.'], ...
        fname, numel(v), n);
end
end

function bt_fwrite(fid, data, fname)
n = numel(data);
if fwrite(fid, data, 'double') ~= n
    error('block_thomas:oocDiskFull', ...
        ['out-of-core write to ''%s'' was truncated: the disk is probably full. ' ...
         'The block-Thomas factor is streamed to disk and needs several GB per ' ...
         'partition at this resolution. Free disk space, lower the resolution, or ' ...
         'use lu_mode=''lapack_band''.'], fname);
end
end
