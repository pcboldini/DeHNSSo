function [Mode,StabRes] = inflow(Mode,Stab,StabGrid,StabRes,Opt)
%INFLOW  Set inflow boundary conditions for all active modes.
%
%   [Mode, StabRes] = inflow(Mode, Stab, StabGrid, StabRes)
%
%   For each mode with nonzero initial amplitude (Stab.A0), computes the
%   inflow condition using ILST (local stability eigenvalue problem) or
%   applies wall/loaded BCs. Sets phase reference (wall pressure or max u).
%   Populates StabRes.u/v/w/p/phi at station 1 and updates Mode.RunJ.

nonzeromodeswall = [];
nonzeromodesic   = [];

% Modes activated by wall/freestream actuation (Stab.bcw / Stab.bct)
nf_bc = min(size(Stab.bcw,3), Mode.nf);
for j = 1:nf_bc
    if any(Stab.bcw(:,:,j),'all') || any(Stab.bct(:,:,j),'all')
        nonzeromodeswall(end+1) = j; %#ok<AGROW>
    end
end

% Initialize solution vector phi = [u v w p]', Amplitude A and streamwise wavenumber alpha 
StabRes.phi    = zeros(Mode.nf,4*StabGrid.ny,StabGrid.nx);
StabRes.alpha  = zeros(Mode.nf,StabGrid.nx);
StabRes.A      = zeros(Mode.nf,StabGrid.nx);
StabRes.u      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.v      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.w      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);
StabRes.p      = zeros(Mode.nf,StabGrid.ny,StabGrid.nx);

if ~exist('Opt','var') || ~isfield(Opt,'bc_top'); bc_top = 'dirichlet'; else; bc_top = Opt.bc_top; end

if Stab.IC=="ILST"
    [~,nonzeromodesic]= find(Stab.A0>0);

    for j = flip(nonzeromodesic) %Run initial condition for all nonzero modes
    if j >= round(Mode.nf/2)% Mode is physical
        [StabRes]=ILST(j,1,StabGrid,StabRes,bc_top); 
    
        % Output from IC_HNS is normalized by max(u'), superimpose amplitude
        StabRes.u(j,:,1) = Stab.A0(j).*StabRes.u(j,:,1);
        StabRes.v(j,:,1) = Stab.A0(j).*StabRes.v(j,:,1);
        StabRes.w(j,:,1) = Stab.A0(j).*StabRes.w(j,:,1);
        StabRes.p(j,:,1) = Stab.A0(j).*StabRes.p(j,:,1);
        StabRes.phi(j,:,1) = [StabRes.u(j,:,1) StabRes.v(j,:,1) StabRes.w(j,:,1) StabRes.p(j,:,1)];
        switch lower(Stab.phaseRef)  % e.g. 'pwall' or 'umax'
        case 'pwall'
            theta0 = angle(StabRes.p(j,end,1));  % wall pressure
        case 'umax'
            [~,iy] = max(abs(StabRes.u(j,:,1))); % max |u'|
            theta0 = angle(StabRes.u(j,iy,1));
        otherwise
            error('Unknown phaseRef option');
        end
        % reference the phase of the inflow perturbation
        StabRes.phi(j,:,1) = StabRes.phi(j,:,1) * exp(-1i*theta0);
    
    else % Apply symmetry for conjugate mode
        j2 = round(Mode.nf/2)-(j-round(Mode.nf/2));% Defines conjugate mode
        
        StabRes.u(j,:,1) = conj(StabRes.u(j2,:,1));
        StabRes.v(j,:,1) = conj(StabRes.v(j2,:,1));
        StabRes.w(j,:,1) = conj(StabRes.w(j2,:,1));
        StabRes.p(j,:,1) = conj(StabRes.p(j2,:,1));
        StabRes.alpha(j,1) = -conj(StabRes.alpha(j,1)); % growth rates should be equal between conjugate modes so -conj
        StabRes.phi(j,:,1) = conj(StabRes.phi(j2,:,1));
    end                                                                            
    end 
    
    if max(abs(Stab.bcw(2:end,:)))==0
    fprintf('Initialising with ILST solution at inflow. \n')
    else
    fprintf('Initialising with ILST solution at inflow and non-homogeneous b.c. at wall. \n')
    end

else

fprintf('Error: no valid inflow condition. \n')
return

end

% nonzero modes are those initiated at inflow + those forced at wall
nonzeromodes = unique([nonzeromodesic, nonzeromodeswall]);
Mode.RunJ = nonzeromodes;            % Used in loops to assess only active modes
Mode.Aactive(Mode.RunJ) = 1;      % Set active modes to introduced modes

% Define phiIC
StabRes.phiIC = squeeze(StabRes.phi(:,:,1));

end

