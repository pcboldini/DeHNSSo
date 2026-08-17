%% Compile banded-LU MEX files (link against MATLAB's LAPACK).
%   - dgbsv_solve        : factor+solve in one call (*gbsv)
%   - dgbsv_factor       : factor only (*gbtrf) -> cacheable band + pivots
%   - dgbsv_solve_cached : solve with a cached factor (*gbtrs)
%   Run once from any working directory.

here = fileparts(mfilename('fullpath'));
srcs = {'dgbsv_solve.c', 'dgbsv_factor.c', 'dgbsv_solve_cached.c'};
for i = 1:numel(srcs)
    src = fullfile(here, srcs{i});
    fprintf('Compiling %s ...\n', src);
    mex('-R2018a', src, '-lmwlapack', '-outdir', here);
end
fprintf('Done. Built into %s\n', here);
