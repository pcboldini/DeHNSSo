function StabRes = postprocess(StabGrid, Init, Mode, Opt, StabRes)
%POSTPROCESS  Normalize shape functions, compute alpha, delete conjugates.

iu = Init.iu;

%% Compute du/dxi
for j = Mode.RunJ
    dux(j,:,:) = -FD1d4o(squeeze(StabRes.u(j,:,:)), StabGrid.dxi) .* StabGrid.xix(end,:);
end

%% Normalize shape functions and compute amplitudes
yinteg = linspace(0, max(StabGrid.eta1D), 4000);

for j = Mode.RunJ
    for i = 1:StabGrid.nx
        [~, k] = max(abs(interp1(StabGrid.etaun, StabRes.u(j,:,i), yinteg, 'spline')));
        dumax = interp1(StabGrid.eta1D, dux(j,:,i), yinteg(k), 'spline');
        umax  = interp1(StabGrid.eta1D, StabRes.u(j,:,i), yinteg(k), 'spline');

        % Alpha a posteriori
        StabRes.alpha(j,i) = (dumax / umax) / iu;

        % Normalize by local u_max
        StabRes.u(j,:,i) = StabRes.u(j,:,i) ./ abs(umax);
        StabRes.v(j,:,i) = StabRes.v(j,:,i) ./ abs(umax);
        StabRes.w(j,:,i) = StabRes.w(j,:,i) ./ abs(umax);
        StabRes.p(j,:,i) = StabRes.p(j,:,i) ./ abs(umax);

        % Amplitude
        if j == round(Mode.nf/2)
            StabRes.A(j,i) = abs(umax);
        else
            StabRes.A(j,i) = 2 * abs(umax);
        end
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
