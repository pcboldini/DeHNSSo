%% Trim the DeHNSSo-v1 step base flow to the step sub-domain.
% Source: heritage Data_Files/BF_SweptWing_Step.mat (232 MB, x = 0.08..517 in lref).
% The v1 paper Fig. 13 window is x = 129.5..215.1 (step at 177.6); keep 125..220.
% Writes bf_cfi_step.mat beside this file (not tracked).
src = '/Users/pcboldini/Desktop/postdoc/codes/heritage/DeHNSSo_original/Data_Files/BF_SweptWing_Step.mat';
s = load(src); BF = s.BF;
keep = BF.X(:).' >= 125 & BF.X(:).' <= 220;
for f = {'X','U','V','W','dxU','dxV','dxW','dyU','dyV','dyW'}
    A = BF.(f{1});
    if isvector(A), BF.(f{1}) = A(keep); else, BF.(f{1}) = A(:, keep); end
end
if size(BF.Y,2) == numel(keep), BF.Y = BF.Y(:, keep); end
fprintf('trimmed: X %.2f..%.2f (n=%d), size(U)=%s, lref=%.4e Re=%.3f\n', min(BF.X), max(BF.X), numel(BF.X), mat2str(size(BF.U)), BF.lref, BF.Re);
save(fullfile(fileparts(mfilename('fullpath')), 'bf_cfi_step.mat'), 'BF', '-v7.3');
