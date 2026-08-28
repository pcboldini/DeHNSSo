function ML = assemble_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
%ASSEMBLE_ML  Per-mode harmonic operator for one (omega, beta) pair.
%
%   ML = ASSEMBLE_ML(omega, beta, MLA1, MLA2, MLA3, MLA4, MLA5, MFD_diff)
%
%   MLA1     — omega- and beta-independent block
%   MLA2     — linear in omega
%   MLA3     — linear in beta
%   MLA4     — quadratic in beta
%   MLA5     — i*beta*d_x cross-viscous block (compressible only; zero sparse
%              for the incompressible operator)
%   MFD_diff — restores no-slip on the (0,0) mean-flow-distortion mode
%
%   Called by solve_modes to build the operator each mode is factorized from,
%   and by backward_error to rebuild the same operator for the residual check.

if omega == 0 && beta == 0
    ML = MLA1 + MFD_diff;
elseif omega == 0
    ML = MLA1 + MLA3*beta + MLA4*beta^2 + MLA5*beta;
elseif beta == 0
    ML = MLA1 + MLA2*omega;
else
    ML = MLA1 + MLA2*omega + MLA3*beta + MLA4*beta^2 + MLA5*beta;
end
end
