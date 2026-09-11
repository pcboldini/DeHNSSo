%% Description
% Shrinks the reduced step base flow to a size that can ship as a gridgen
% benchmark, and writes baseflow/output/benchmark/bf_cfi_step.mat.
%
% Input is bf_cfi_step_dns.mat, itself a pure subsample of the DNS written
% by make_bf_from_dns.m, so the result is still made only of DNS nodes and no
% value is ever interpolated. Working from it avoids fetching the 1.8 GB source
% back from lsl; it keeps full DNS resolution below y = 2.5 and within 6 of the
% face, which is finer than every target below.
%
% Coarsening was tuned against the stability grid the benchmark caller builds
% (1500 x 80, y_i = 1). Interpolation error is concentrated in one place: the
% shear layer springing off the step edge, within about one lref downstream of
% the face at y just above the step top. Full DNS resolution is therefore kept
% for y <= 2 and within 2 of the face, and the grid coarsens steeply outside,
% where a first attempt at 0.03 and 0.027 left errors below 1e-3 anyway.
% Near the wall the stability grid is far finer than any base flow (dy = 4.1e-4)
% because of the Malik clustering, but the profile is nearly linear there.

clear all; close all; clc
here = fileparts(mfilename('fullpath'));
root = fullfile(here, '..', '..', '..');

L = load(fullfile(here, 'bf_cfi_step_dns.mat'));  B = L.BF;
xf = B.StepX;  hs = B.StepH;
fprintf('input %d x %d (%.1f MB)\n', numel(B.X), numel(B.Y), ...
        dir(fullfile(here,'bf_cfi_step_dns.mat')).bytes/1e6);

%% Target spacings
x = B.X(:).';   y = B.Y(:);            % y descending
dxt = repmat(0.25, numel(x), 1);
dxt(abs(x - xf) <= 20) = 0.06;
dxt(abs(x - xf) <=  2) = 0;        % 0 keeps every DNS column
keepx = greedy_keep(x(:), dxt);
ifac  = find(x >= xf - 1e-9, 1);
keepx([ifac-1, ifac, 1, numel(x)]) = true;

ya  = flipud(y);                        % ascending for the rule
dyt = repmat(1.00, numel(ya), 1);
dyt(ya <= 10) = 0.25;
dyt(ya <=  6) = 0.08;
dyt(ya <=  2) = 0;                 % 0 keeps every DNS row
keepya = greedy_keep(ya, dyt);
[~, jh] = min(abs(ya - hs));  keepya(jh) = true;      % the node on the step top
keepya([1 numel(ya)]) = true;
keepy = flipud(keepya);

%% Package
BF.X = x(keepx);  BF.Y = y(keepy);
for f = {'U','V','W','dxU','dxV','dxW','dyU','dyV','dyW'}
    BF.(f{1}) = B.(f{1})(keepy, keepx);
end
BF.lref = B.lref;  BF.Uref = B.Uref;  BF.Re = B.Re;  BF.nu = B.nu;
BF.StepX = B.StepX;  BF.StepH = B.StepH;
BF.source = 'dataExtCrossD.mat (lsl:~/postdoc/dns_data), spanwise mean, DNS nodes only';

out = fullfile(root, 'baseflow', 'output', 'benchmark', 'bf_cfi_step.mat');
save(out, 'BF', '-v7.3');
fprintf('output %d x %d (%.1f MB, %.1fx fewer points)\n', numel(BF.X), numel(BF.Y), ...
        dir(out).bytes/1e6, numel(B.X)*numel(B.Y)/(numel(BF.X)*numel(BF.Y)));
dxb = abs(diff(BF.X));  dyb = abs(diff(BF.Y));
fprintf('  dx %.4f .. %.4f,  dy %.4f .. %.4f\n', min(dxb), max(dxb), min(dyb), max(dyb));
fprintf('  node on the step top kept: %d,  solid still exactly zero: %d\n', ...
        any(abs(BF.Y - hs) < 1e-9), max(abs(BF.U(BF.Y > 0 & BF.Y <= hs, BF.X >= xf)), [], 'all') == 0);

function k = greedy_keep(s, target)
% Keep the first point, then every point at least target(i) beyond the last kept.
n = numel(s);  k = false(n,1);  k(1) = true;  last = s(1);
for i = 2:n
    if target(i) <= 0 || abs(s(i) - last) >= target(i)
        k(i) = true;  last = s(i);
    end
end
end
