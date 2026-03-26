function print_config(StabGrid, Stab, Opt)
%PRINT_CONFIG  Print run configuration and memory estimates.

fprintf('  --------------------------------------------\n');
fprintf('  Grid\n');
fprintf('    nx = %d,  ny = %d,  DOF = %d\n', StabGrid.nx, StabGrid.ny, 4*StabGrid.ny*StabGrid.nx);
fprintf('    Re = %.2f\n', StabGrid.Re);
fprintf('    H = %.4f,  y_i = %.4f\n', StabGrid.H, StabGrid.y_i);
fprintf('    x = [%.2f, %.2f],  dxi = %.4f\n', StabGrid.x(1), StabGrid.x(end), StabGrid.dxi);
fprintf('  Spectral content\n');
fprintf('    M = %d,  N = %d  (%d modes)\n', Stab.M, Stab.N, (2*Stab.M+1)*(2*Stab.N+1));
fprintf('    omega_0 = %.6f,  beta_0 = %.6f\n', Stab.omega_0, Stab.beta_0);
fprintf('  Inflow\n');
fprintf('    IC method: %s,  active modes: %d\n', Stab.IC, nnz(Stab.A0));
fprintf('  Solver\n');
fprintf('    Buffer: %s from %d%%,  Conv: %.1e\n', Opt.buffer, Opt.xb, Opt.Conv);
fprintf('    LU mode: %s\n', Opt.lu_mode);
if isfield(Opt,'gpu') && strcmpi(Opt.gpu,'on')
    try
        g = gpuDevice;
        fprintf('    GPU: %s (%.1f GB VRAM)\n', g.Name, g.TotalMemory/1e9);
    catch
        fprintf('    GPU: enabled (device info unavailable)\n');
    end
end

% Memory estimates
nx = StabGrid.nx;  ny = StabGrid.ny;  bsize = 4*ny;
N_dof = bsize * nx;
D1_is_sparse = isfield(StabGrid,'D1') && issparse(StabGrid.D1) && (nnz(StabGrid.D1) < 0.2*ny^2);
if D1_is_sparse
    lhs_nnz_est = (20 + 2) * N_dof;
    eta_lbl = 'fd';
else
    lhs_nnz_est = (54 * (ny/40) + 2) * N_dof;
    eta_lbl = 'cheb';
end
lhs_mem_est = lhs_nnz_est * 20 / 1e9;

% LU fill-in is bandwidth-governed: L and U each fill a band of width b_half.
% b_half = 2*bsize for the 5-point xi stencil (spans ±2 stations).
% AMD reordering reduces fill by ~25% in practice (factor 0.75).
% 20 bytes/nnz (complex double: 16 bytes + 4 bytes index overhead).
b_half      = 2 * bsize;
lu_mem_one  = 2 * N_dof * b_half * 20 * 0.75 / 1e9;
n_unique    = Stab.M + 1;
if ~isfield(Opt, 'lu_max_cache'); Opt.lu_max_cache = 3; end

% parfor workers each hold a copy of the broadcast LU objects
% nw: number of parallel workers (default pool size, capped at n_unique modes)
if strcmpi(Opt.parfor,'on') && license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'))
    pool = gcp('nocreate');
    if isempty(pool); nw = 0; else; nw = min(pool.NumWorkers, n_unique); end
else
    nw = 0;
end
parfor_overhead = nw * n_unique * lu_mem_one;  % each worker gets all LU objects

fprintf('    Memory estimates (%s):\n', eta_lbl);
fprintf('      LHS:                   ~%.1f GB\n', lhs_mem_est);
fprintf('      LU per factorisation:  ~%.1f GB\n', lu_mem_one);
fprintf('      LU full   (cache %d):  ~%.1f GB\n', n_unique, n_unique * lu_mem_one);
fprintf('      LU auto   (cache %d):  ~%.1f GB\n', Opt.lu_max_cache, Opt.lu_max_cache * lu_mem_one);
fprintf('      LU single (cache 1):   ~%.1f GB\n', lu_mem_one);
fprintf('      Total peak (single):   ~%.1f GB\n', lu_mem_one + lhs_mem_est);
if nw > 0
    fprintf('      parfor workers (%d):   ~%.1f GB overhead\n', nw, parfor_overhead);
    fprintf('      Total peak (parfor):   ~%.1f GB\n', n_unique * lu_mem_one + lhs_mem_est + parfor_overhead);
end

lu_mem_total = lu_mem_one;
if strcmp(Opt.lu_mode, 'full'); lu_mem_total = n_unique * lu_mem_one; end
if strcmpi(Opt.parfor,'on'); lu_mem_total = lu_mem_total + parfor_overhead; end
if lu_mem_total > 20
    warning('Estimated LU memory: %.1f GB. Consider lu_mode = ''auto'' or ''none'' if RAM is limited.', lu_mem_total);
end
fprintf('  --------------------------------------------\n');

end
