%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2026 P.C. Boldini, Y. Chen & M. Kotsonis
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
%
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
%
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <https://www.gnu.org/licenses/>.
%
%     Please cite this code and the paper if you have used this for your
%     publication:

% DOI:

%% Description
% Reduces the three-dimensional DNS of the forward-facing-step swept-wing case
% to the two files this test case needs, following the procedure used for the
% curved swept-wing case of Sec. 6.3:
%
%   1. base flow  = spanwise mean of the DNS velocity field and its derivatives,
%                   restricted to the stability window and coarsened to the
%                   resolution the stability grid can use;
%   2. reference  = the spanwise Fourier mode amplitudes A_n(x) precomputed in
%                   the DNS files as AMaxU, for the stepped and the clean plate.
%
% Source files (read-only, not in the repository, and not on this machine: the
% laptop volume has no room for them. Both live on lsl in ~/postdoc/dns_data/,
% so point FD/FA there or copy them back to re-run this. The two outputs written
% below are what the case itself needs):
%   dataExtCrossD.mat   forward-facing step
%   dataExtCrossA.mat   clean plate, same actuator and spectral content
% Both hold g.u, g.v, g.w and their gradients on 6760 x 72 x 576 points, which
% is 2.8e8 points per field and 31 GB if read in full. Only the window
% 125 <= x <= 220 is read, one field at a time, in streamwise chunks.
%
% The step geometry is imposed from the DNS metadata g.xLoc and g.h rather than
% inferred: the face column carries a nonzero wall-normal velocity in the DNS,
% so a detector keyed on all three components places the face half a cell
% downstream of the reported location.
%
% Writes bf_cfi_step_dns.mat and ref_dns_step_modes.mat in this folder.

clear all; close all; clc

here = fileparts(mfilename('fullpath'));
FD   = '/Users/pcboldini/Downloads/dataExtCrossD.mat';   % step
FA   = '/Users/pcboldini/Downloads/dataExtCrossA.mat';   % clean
if ~isfile(FD) || ~isfile(FA)
    error('make_bf_from_dns:noSource', ...
        ['The DNS source files are not present. Fetch dataExtCrossD.mat and ' ...
         'dataExtCrossA.mat (about 2 GB each) and point FD/FA at them. The reduced ' ...
         'outputs bf_cfi_step_dns.mat and ref_dns_step_modes.mat are already ' ...
         'in this folder, so the case runs without them.']);
end

xwin = [125, 220];        % base-flow window, wider than the stability domain
dx_near = 0.0135;         % target spacing within x_face +- 6   (DNS spacing)
dx_mid  = 0.040;          % target spacing within x_face +- 25
dx_far  = 0.100;          % target spacing elsewhere
y_fine  = 2.5;            % keep every DNS point below this height
dy_mid  = 0.060;          % target spacing for y_fine < y <= 8
dy_far  = 0.300;          % target spacing above y = 8

%% Geometry and coordinates
Xun  = h5read(FD, '/g/Xun');   Xun = Xun(:);        % ascending
Yun  = h5read(FD, '/g/Yun');   Yun = Yun(:);        % descending, 25.86 -> 0.0068
lref = h5read(FD, '/g/lref');  uref = h5read(FD, '/g/uref');
Reg  = h5read(FD, '/g/Reg');
xLoc = h5read(FD, '/g/xLoc');  hStep = h5read(FD, '/g/h');
nz   = numel(h5read(FD, '/g/Zun'));

iw  = find(Xun >= xwin(1) & Xun <= xwin(2));
i0  = iw(1);  nxw = numel(iw);  ny = numel(Yun);
fprintf('DNS %d x %d x %d; window x = %.2f..%.2f -> %d stations\n', ...
        numel(Xun), nz, ny, Xun(i0), Xun(iw(end)), nxw);
fprintf('step from metadata: x_face = %.6f, h = %.6f\n', xLoc, hStep);

%% Spanwise mean of every field, read in streamwise chunks
src = {'u','v','w','ux','vx','wx','uy','vy','wy'};
dst = {'U','V','W','dxU','dxV','dxW','dyU','dyV','dyW'};
chunk = 250;
for f = 1:numel(src)
    A = zeros(ny, nxw);
    for c = 1:chunk:nxw
        n = min(chunk, nxw - c + 1);
        blk = h5read(FD, ['/g/' src{f}], [1 1 i0+c-1], [ny nz n]);
        A(:, c:c+n-1) = squeeze(mean(blk, 2));
    end
    BFfull.(dst{f}) = A;
    fprintf('  spanwise mean %-3s -> %s done\n', src{f}, dst{f});
end
clear blk A

%% Impose the solid region from the DNS metadata
% g.xLoc falls exactly on a streamwise node and g.h within 5e-9 of a wall-normal
% node, so both are snapped to those nodes: the DNS itself carries no-slip on the
% surface row, which is the top of the solid.
Xw = Xun(iw);
[~, jtop]  = min(abs(Yun - hStep));   hnode = Yun(jtop);
[~, iface] = min(abs(Xw  - xLoc));    xnode = Xw(iface);
fprintf('snapped to nodes: face x = %.9f (%.1e off), top y = %.9f (%.1e off)\n', ...
        xnode, xnode - xLoc, hnode, hnode - hStep);
solid = logical((Yun(:) > 0 & Yun(:) <= hnode) * (Xw(:).' >= xnode));
hStep = hnode;   xLoc = xnode;
fprintf('solid cells zeroed: %d of %d (%.2f%%)\n', nnz(solid), numel(solid), 100*nnz(solid)/numel(solid));
for f = 1:numel(dst); BFfull.(dst{f})(solid) = 0; end

%% Point reduction
dxt = repmat(dx_far, nxw, 1);
dxt(abs(Xw - xLoc) <= 25) = dx_mid;
dxt(abs(Xw - xLoc) <=  6) = dx_near;
keepx = greedy_keep(Xw, dxt);
ifac  = iface;                                      % first solid column
keepx([ifac-1, ifac, 1, nxw]) = true;               % keep the face bracket and the ends

Ya   = flipud(Yun);                                  % ascending for the rule
dyt  = repmat(dy_far, ny, 1);
dyt(Ya <= 8)      = dy_mid;
dyt(Ya <= y_fine) = 0;                               % 0 -> keep every point
keepya = greedy_keep(Ya, dyt);
[~, jh] = min(abs(Ya - hStep));  keepya(jh) = true;  % keep the node on the step top
keepya([1 ny]) = true;
keepy = flipud(keepya);                              % back to the stored order

nxr = nnz(keepx);  nyr = nnz(keepy);
fprintf('reduction: %d x %d -> %d x %d  (%.2fM -> %.2fM points, %.1fx)\n', ...
        nxw, ny, nxr, nyr, nxw*ny/1e6, nxr*nyr/1e6, (nxw*ny)/(nxr*nyr));

BF.X = Xw(keepx).';  BF.Y = Yun(keepy);
for f = 1:numel(dst); BF.(dst{f}) = BFfull.(dst{f})(keepy, keepx); end
BF.lref = lref;  BF.Uref = uref;  BF.Re = Reg;  BF.nu = uref*lref/Reg;
BF.StepX = xLoc; BF.StepH = hStep;

%% Check against the heritage base flow, which is the same spanwise mean unreduced
her = fullfile(here, 'bf_cfi_step.mat');
if isfile(her)
    H = load(her);  Hb = H.BF;
    [~, ix] = ismembertol(BF.X(:), Hb.X(:), 1e-9, 'DataScale', 1);
    [~, iy] = ismembertol(BF.Y(:), Hb.Y(:), 1e-9, 'DataScale', 1);
    ok = ix > 0;  jy = iy > 0;
    for f = {'U','V','W'}
        a = BF.(f{1})(jy, ok);  b = Hb.(f{1})(iy(jy), ix(ok));
        d = abs(a - b);
        fprintf('  vs heritage %s: max|diff| outside solid = %.3e, inside = %.3e (max|%s| = %.4f)\n', ...
                f{1}, max(d(~solid(jy, ok))), max([d(solid(jy, ok)); 0]), f{1}, max(abs(b(:))));
    end
end

save(fullfile(here, 'bf_cfi_step_dns.mat'), 'BF', '-v7.3');
d = dir(fullfile(here, 'bf_cfi_step_dns.mat'));
fprintf('wrote bf_cfi_step_dns.mat (%.1f MB)\n', d.bytes/1e6);

%% Reference mode amplitudes, both configurations
ref = struct();
for F = {{'step', FD}, {'clean', FA}}
    tag = F{1}{1};  file = F{1}{2};
    X = h5read(file, '/g/Xun'); X = X(:);
    A = h5read(file, '/g/AMaxU');            % (6 modes, nx)
    if size(A,2) ~= numel(X); A = A.'; end
    ref.(tag).x    = X;
    ref.(tag).A    = A;                      % row 1 = MFD (0,0), rows 2..6 = (0,1)..(0,5)
    ref.(tag).modes = 0:5;
end
ref.lref = lref; ref.uref = uref; ref.Re = Reg;
ref.x_paper = @(x) 3.604*x + 220.6;
ref.note = 'A(1,:) = MFD (0,0); A(n+1,:) = harmonic (0,n), normalised by uref';
save(fullfile(here, 'ref_dns_step_modes.mat'), '-struct', 'ref', '-v7');
fprintf('wrote ref_dns_step_modes.mat\n');
for tag = {'step','clean'}
    A = ref.(tag{1}).A; X = ref.(tag{1}).x;
    q = @(xt) A(:, find(abs(X-xt) == min(abs(X-xt)), 1));
    fprintf('%-6s at x=150: (0,0) %.3e (0,1) %.3e | at x=%.1f: (0,0) %.3e (0,1) %.3e\n', ...
            tag{1}, sub(q(150),1), sub(q(150),2), xLoc, sub(q(xLoc),1), sub(q(xLoc),2));
end

function k = greedy_keep(s, target)
% Keep the first point, then every point at least target(i) beyond the last kept.
% target = 0 keeps every point.
n = numel(s);  k = false(n,1);  k(1) = true;  last = s(1);
for i = 2:n
    if target(i) <= 0 || abs(s(i) - last) >= target(i)
        k(i) = true;  last = s(i);
    end
end
end

function v = sub(a, i); v = a(i); end
