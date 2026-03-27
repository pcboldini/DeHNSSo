function StabRes = postprocess(StabGrid, Init, Mode, Opt, StabRes)
%POSTPROCESS  Normalize shape functions, compute alpha, delete conjugates.

iu = Init.iu;

%% Compute du/dxi
for j = Mode.RunJ
    dux(j,:,:) = -FD1d4o(squeeze(StabRes.u(j,:,:)), StabGrid.dxi) .* StabGrid.xix(end,:);
end

%% Normalize shape functions and compute amplitudes
% Interpolate all stations at once using matrix interp1 — one spline per mode
% instead of 2*nx scalar spline calls per mode.
yinteg = linspace(0, max(StabGrid.eta1D), 4000);

for j = Mode.RunJ
    % Interpolate u and dux over all stations simultaneously [4000 x nx]
    U_fine   = interp1(StabGrid.etaun, squeeze(StabRes.u(j,:,:)), yinteg, 'spline');
    DUX_fine = interp1(StabGrid.eta1D, squeeze(dux(j,:,:)),       yinteg, 'spline');

    % Find wall-normal location of max |u| per station
    [~, k_all] = max(abs(U_fine), [], 1);  % 1 x nx

    % Extract umax and dumax at the per-station maximum location
    idx      = sub2ind(size(U_fine), k_all, 1:StabGrid.nx);
    umax_vec  = U_fine(idx);
    dumax_vec = DUX_fine(idx);

    % Alpha a posteriori (vectorised)
    StabRes.alpha(j,:) = (dumax_vec ./ umax_vec) / iu;

    % Normalise shape functions and compute amplitudes (vectorised over stations)
    absU = abs(umax_vec);                           % 1 x nx
    StabRes.u(j,:,:) = squeeze(StabRes.u(j,:,:)) ./ absU;
    StabRes.v(j,:,:) = squeeze(StabRes.v(j,:,:)) ./ absU;
    StabRes.w(j,:,:) = squeeze(StabRes.w(j,:,:)) ./ absU;
    StabRes.p(j,:,:) = squeeze(StabRes.p(j,:,:)) ./ absU;

    if j == round(Mode.nf/2)
        StabRes.A(j,:) = absU;
    else
        StabRes.A(j,:) = 2 * absU;
    end
end

%% Delete conjugates
half = round(Mode.nf/2);
StabRes.omegavec(:, 1:half-1) = [];
StabRes.betavec(:, 1:half-1)  = [];
StabRes.A(1:half-1, :)        = [];
StabRes.alpha(1:half-1, :)    = [];
StabRes.u(1:half-1, :, :)     = [];
StabRes.v(1:half-1, :, :)     = [];
StabRes.w(1:half-1, :, :)     = [];
StabRes.p(1:half-1, :, :)     = [];

if Opt.Sweep
    StabRes.Asweep(1:half-1, :, :)       = [];
    StabRes.alphasweep(1:half-1, :, :)   = [];
    StabRes.usweep(1:half-1, :, :, :)    = [];
    StabRes.vsweep(1:half-1, :, :, :)    = [];
    StabRes.wsweep(1:half-1, :, :, :)    = [];
    StabRes.psweep(1:half-1, :, :, :)    = [];
end

%% Remove working arrays
StabRes = rmfield(StabRes, {'q', 'phi', 'R'});

end
