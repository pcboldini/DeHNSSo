function [Opt,StabRes] = AVAL(Grid,Opt,Mode,StabRes)
%AVAL  Check the peak disturbance amplitude and damp the inflow if it is too large.
%
%   [Opt, StabRes] = AVAL(Grid, Opt, Mode, StabRes)
%
%   Finds the peak velocity amplitude of the fundamental mode and stores it in
%   Opt.AVAL. If that peak exceeds the ceiling Opt.AMAX, the amplitude factor
%   Opt.AF is rescaled so the peak comes back down to Opt.AMAX, and the whole
%   solution vector is scaled by the new Opt.AF.

%% Reshape the solution vector into the field array phi
    ny_   = Grid.ny;
    nx_   = Grid.nx;
    nvar  = size(StabRes.phi, 2) / ny_;
    n_vel = 3;                              % u, v, w slots
    for j = Mode.RunJ
        StabRes.phi(j,:,:) = reshape(StabRes.q(j,:), nvar*ny_, nx_);
    end

%% Peak velocity amplitude of the fundamental mode
    [Opt.AVAL,~] = max(max(abs(squeeze(StabRes.phi(round(Mode.nf/2)+1,1:n_vel*ny_,:)))));

%% Damp the inflow amplitude if the peak exceeds the ceiling AMAX
    if Opt.AVAL >= Opt.AMAX
        Opt.AF = Opt.AMAX/Opt.AVAL;
    end

    StabRes.q = Opt.AF*StabRes.q;
end
