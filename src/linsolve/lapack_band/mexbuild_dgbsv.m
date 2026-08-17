%% Compile the banded LU (LU factorization) MEX (MATLAB executable, compiled
%   C) files, linking against MATLAB's own LAPACK.
%   - dgbsv_solve        : factor and solve in one call (*gbsv)
%   - dgbsv_factor       : factor only (*gbtrf) -> a cacheable band matrix + pivots
%   - dgbsv_solve_cached : solve reusing a cached factor, no re-factoring (*gbtrs)
%   Run this once, from any working directory.

here = fileparts(mfilename('fullpath'));
srcs = {'dgbsv_solve.c', 'dgbsv_factor.c', 'dgbsv_solve_cached.c'};
for i = 1:numel(srcs)
    src = fullfile(here, srcs{i});
    fprintf('Compiling %s ...\n', src);
    mex('-R2018a', src, '-lmwlapack', '-outdir', here);
end
fprintf('Done. Built into %s\n', here);
