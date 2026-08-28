function B = bc_rows(ny, is_compr)
%BC_ROWS  Row indices, within one station block, that carry the boundary conditions.
%
%   B = bc_rows(ny, is_compr)
%
%   Single source of truth for which row of a station block each boundary-condition
%   (BC) channel k acts on. LHS/LHS_compressible impose the BC type on these rows and
%   RHS writes the inhomogeneous values (Stab.bcw at the wall, Stab.bct at the
%   freestream) into the same rows, so the two must agree channel for channel.
%
%   Grid convention: index 1 = freestream (y=H), index ny = wall (y=0). The state is
%   stacked (u,v,w,p) incompressible and (u,v,w,rho,T) compressible, each ny long.
%
%   Fields:
%     B.top        freestream row per channel
%     B.top_mfd    subset of B.top also constrained for the mean-flow distortion (MFD)
%     B.wall       wall row per channel
%     B.top_names  variable each channel constrains at the freestream
%     B.wall_names variable each channel constrains at the wall
%
%   NOTE the channel maps are NOT the same at the two boundaries in the compressible
%   case, because the two boundaries constrain different variables:
%     freestream  u, w, rho, T   (v is left free: it carries the displacement)
%     wall        u, v, w, T     (rho is left free: it follows from the equations)
%   Incompressibly both are (u,v,w) and p is free at the freestream.

if is_compr
    B.top        = [1, 2*ny+1, 3*ny+1, 4*ny+1];   % u, w, rho, T
    B.top_mfd    = B.top;
    B.wall       = [ny, 2*ny, 3*ny, 5*ny];        % u, v, w, T
    B.top_names  = {'u','w','rho','T'};
    B.wall_names = {'u','v','w','T'};
else
    B.top        = [1, ny+1, 2*ny+1];             % u, v, w
    B.top_mfd    = [1, 2*ny+1];                   % MFD v is free (displacement velocity)
    B.wall       = [ny, 2*ny, 3*ny];              % u, v, w
    B.top_names  = {'u','v','w'};
    B.wall_names = {'u','v','w'};
end
end
