function [Bound]=setup_boundary(Grid,Opt,Stab)
%SETUP_BOUNDARY  Set up the outflow buffer zone and inhomogeneous boundary conditions (BCs).
%
%   Bound = setup_boundary(Grid, Opt, Stab)
%
%   Builds the buffer region near the outflow, where the solution is smoothly
%   damped or parabolized to prevent reflections. The buffer ramp starts at
%   Opt.xb (percent of the domain) and steepens with Opt.kappa. Also
%   interpolates the inhomogeneous wall/freestream BCs (Stab.bcw, Stab.bct)
%   from wherever they were defined onto the numerical grid and stores the
%   result in Bound.bcw, Bound.bct.
%
%   Output: Bound struct with fields bufc, bufs, bufsp, bufcf, ibuf,
%           bcw (3 x nx x nf), bct (3 x nx x nf).

%% Buffer

switch lower(Opt.buffer)
case 'on'
    % Amplitude-damping buffer: scales the whole finite-difference (FD) block by bufc -> 0
    Bound.ibuf  = round(Opt.xb * Grid.nx / 100);
    Bound.bufc  = ones(1, Grid.nx);
    Bound.bufp  = ones(1, Grid.nx);   % Fi (streamwise diffusion term) is left unchanged in 'on' mode
    Bound.bufc  = buffer_ramp(Bound.bufc,  Bound.ibuf,  Grid.nx, Opt.kappa);

    % NLT forcing buffer
    Bound.ibuff = round(Opt.nltbufxb * Grid.nx / 100);
    Bound.bufcf = ones(1, Grid.nx);
    Bound.bufcf = buffer_ramp(Bound.bufcf, Bound.ibuff, Grid.nx, Opt.kappa);

    fprintf(' Buffer: amplitude damping from %.3g%% to 100%%\n', Opt.xb);

case 'para'
    % Parabolization buffer: ramps Fi (streamwise diffusion) to zero in the buffer,
    % keeping advection and wall-normal diffusion. This makes the outflow behave
    % like a linear parabolized stability equation (LPSE) solve, avoiding reflections.
    Bound.ibuf  = round(Opt.xb * Grid.nx / 100);
    Bound.bufc  = ones(1, Grid.nx);   % no amplitude damping in this mode

    Bound.bufp  = ones(1, Grid.nx);
    Bound.bufp  = buffer_ramp(Bound.bufp,  Bound.ibuf,  Grid.nx, Opt.kappa);

    % NLT forcing buffer (same ramp as bufp)
    Bound.ibuff = round(Opt.nltbufxb * Grid.nx / 100);
    Bound.bufcf = ones(1, Grid.nx);
    Bound.bufcf = buffer_ramp(Bound.bufcf, Bound.ibuff, Grid.nx, Opt.kappa);

    fprintf(' Buffer: parabolisation from %.3g%% to 100%%\n', Opt.xb);

case 'off'
    Bound.ibuf  = Grid.nx + 1;
    Bound.ibuff = Grid.nx + 1;
    Bound.bufc  = ones(1, Grid.nx);
    Bound.bufp  = ones(1, Grid.nx);
    Bound.bufcf = ones(1, Grid.nx);
    fprintf(' Buffer: off\n');
end

%% Wall masks (this release does not support sharp wall features)
Bound.bufs  = ones(Grid.ny, Grid.nx);
Bound.bufsp = ones(Grid.ny, Grid.nx);
Bound.stepi = 0;

%% Inhomogeneous BCs: interpolate Stab.bcw / Stab.bct onto numerical grid
is_compr = isfield(Opt,'compressible') && strcmpi(Opt.compressible,'on');
n_bc     = 3 + double(is_compr);   % u,v,w [,T]
nf = size(Stab.bcw, 3);
Bound.bcw = zeros(n_bc, Grid.nx, nf);
Bound.bct = zeros(n_bc, Grid.nx, nf);
has_bcw = any(Stab.bcw(:));
has_bct = any(Stab.bct(:));
if has_bcw || has_bct
    for j = 1:nf
        for k = 1:n_bc
            if has_bcw && any(Stab.bcw(k,:,j))
                Bound.bcw(k,:,j) = interp1(Stab.bcx, Stab.bcw(k,:,j), Grid.xun, 'cubic', 0);
            end
            if has_bct && any(Stab.bct(k,:,j))
                Bound.bct(k,:,j) = interp1(Stab.bcx, Stab.bct(k,:,j), Grid.xun, 'cubic', 0);
            end
        end
    end
    if has_bcw; fprintf(' Inhomogeneous wall BC active (%d modes).\n', nnz(squeeze(any(any(Stab.bcw,1),2)))); end
    if has_bct; fprintf(' Inhomogeneous freestream BC active (%d modes).\n', nnz(squeeze(any(any(Stab.bct,1),2)))); end
end

end

%% Local helpers

function buf = buffer_ramp(buf, ib, nx, kappa)
%BUFFER_RAMP  Apply a tanh ramp over buf(ib:nx), normalized to 1 at ib.
i_vec        = ib:nx;
span         = max(nx - ib, 1);
buf(i_vec)   = 0.5*(1 + tanh(kappa*(1 - 2*(i_vec - ib)/span)));
buf(i_vec)   = buf(i_vec) ./ buf(ib);
end

