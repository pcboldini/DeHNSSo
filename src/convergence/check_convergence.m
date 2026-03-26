function [Opt, Stab, StabRes] = check_convergence(StabGrid, Bound, Mode, Opt, Stab, StabRes)
%CHECK_CONVERGENCE  Evaluate convergence, ramp amplitude, extrapolate forcing.

%% Convergence measure
for j = Mode.RunJ
    Opt.dal(j) = sum(abs(StabRes.A(j,1:Bound.ibuf) - StabRes.Aold(j,1:Bound.ibuf))) ...
               / max(StabRes.A(j,1:Bound.ibuf)) / StabGrid.nx;
end

% Linear simulation: force convergence after first iteration
if length(Mode.RunJ) == 1
    Opt.dal = zeros(size(Opt.dal));
end

% Divergence check
if Opt.AVAL >= 2 && Opt.k1 ~= 0
    Opt.dal = zeros(size(Opt.dal));
end

% Save amplitudes
StabRes.Aold = StabRes.A;

%% Amplitude ramping
if Opt.k1 ~= 0
    if max(abs(Opt.dal)) <= (1-Opt.Sweep)*Opt.Conv*Opt.ConvF + Opt.Sweep*Opt.Conv
        Opt.AF = Opt.AF * Opt.AFg;
        if ge(Opt.AF, 1); Opt.AF = 1; end

        % Initialize forcing extrapolation list
        if ~isfield(Opt, 'AFlist')
            Opt.AFlist = [Opt.AF Opt.AF];
            Stab.flist = Stab.f;
            Stab.flist(:,:,size(Stab.flist,3)+1) = Stab.f;
        end

        % Append
        Opt.AFlist = [Opt.AFlist Opt.AF];
        Stab.flist(:,:,size(Stab.flist,3)+1) = Stab.f;

        % Linearly (in A^2) extrapolate the forcing term
        if Opt.AFlist(end) ~= Opt.AFlist(end-1) && Opt.AFlist(end-1) ~= Opt.AFlist(end-2)
            Stab.f = Stab.flist(:,:,end) + ...
                (Stab.flist(:,:,end) - Stab.flist(:,:,end-1)) / ...
                (Opt.AFlist(end-1)^2 - Opt.AFlist(end-2)^2) * ...
                (Opt.AFlist(end)^2 - Opt.AFlist(end-1)^2);
        end
    end
end

%% Update iteration counter and convergence log
Opt.k1 = Opt.k1 + 1;
Opt.dalsave(:, Opt.k1) = abs(Opt.dal);
Opt.dalsave(Opt.dalsave <= 0) = NaN;
Opt.dalsave(Opt.dalsave == 1) = NaN;

end
