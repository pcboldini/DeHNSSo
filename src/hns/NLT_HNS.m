function [ f ] = NLT_HNS(Grid,RunJ,StabRes,D1,HB)
%% License (GNU GENERAL PUBLIC LICENSE v3)
%                  Delft Harmonic Navier-Stokes Solver
%     Copyright (C) 2023 S.H.J. Westerbeek, S. Hulshoff, H. Schuttelaars
%                          & M. Kotsonis
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
% DOI:

%% Find dimensions
[nf,ny,nx] = size(StabRes.u);
iu = 1i;

% Reset f
f = zeros(nf, 4*ny*nx);

% Metric coefficients
etay = Grid.etay;
etax = Grid.etax;
xiy  = Grid.xiy;
xix  = Grid.xix;
dxi  = Grid.x(2) - Grid.x(1);

%% Calculate derivative fields (batched over modes)

% xi-derivatives: apply FD to each active mode's 2D field
dudxi = zeros(nf, ny, nx);
dvdxi = zeros(nf, ny, nx);
dwdxi = zeros(nf, ny, nx);
for j = RunJ
    dudxi(j,:,:) = -FD1d4o(squeeze(StabRes.u(j,:,:)), dxi);
    dvdxi(j,:,:) = -FD1d4o(squeeze(StabRes.v(j,:,:)), dxi);
    dwdxi(j,:,:) = -FD1d4o(squeeze(StabRes.w(j,:,:)), dxi);
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

%% Calculate forcing (accumulate directly, no intermediate arrays)

xmom = zeros(nf, ny, nx);
ymom = zeros(nf, ny, nx);
zmom = zeros(nf, ny, nx);

for j = RunJ(RunJ >= round(nf/2))

    % Find interacting mode pairs from Harmonic Balancing
    [gharray, jkarray] = find(HB(:,:,j));

    % Keep only active mode pairs
    active = ismember(gharray, RunJ) & ismember(jkarray, RunJ);
    gharray = gharray(active);
    jkarray = jkarray(active);

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

        % Accumulate momentum forcing directly (no temporary arrays stored)
        xmom(j,:,:) = squeeze(xmom(j,:,:)) ...
            - u_gh.*dudxi_jk.*xix - u_gh.*dudeta_jk.*etax ...
            - v_gh.*dudxi_jk.*xiy - v_gh.*dudeta_jk.*etay ...
            - iu*beta .* w_gh.*u_jk;

        ymom(j,:,:) = squeeze(ymom(j,:,:)) ...
            - u_gh.*dvdxi_jk.*xix - u_gh.*dvdeta_jk.*etax ...
            - v_gh.*dvdxi_jk.*xiy - v_gh.*dvdeta_jk.*etay ...
            - iu*beta .* w_gh.*v_jk;

        zmom(j,:,:) = squeeze(zmom(j,:,:)) ...
            - u_gh.*dwdxi_jk.*xix - u_gh.*dwdeta_jk.*etax ...
            - v_gh.*dwdxi_jk.*xiy - v_gh.*dwdeta_jk.*etay ...
            - iu*beta .* w_gh.*w_jk;
    end

    % Reshape momentum into f (vectorised, skip i=1)
    xm = squeeze(xmom(j,:,2:nx));  % [ny x nx-1]
    ym = squeeze(ymom(j,:,2:nx));
    zm = squeeze(zmom(j,:,2:nx));
    zr = zeros(ny, nx-1);
    block = [xm; ym; zm; zr];      % [4*ny x nx-1]
    f(j, 4*ny+1:end) = block(:).';

end

end
