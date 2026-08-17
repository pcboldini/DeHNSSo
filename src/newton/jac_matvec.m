function w = jac_matvec(v_real, A, B, N)
%JAC_MATVEC  Real-arithmetic matrix-vector product (matvec) for the Newton
%   Jacobian of the HNS (harmonic Navier-Stokes) system.
%
%   w = jac_matvec(v_real, A, B, N)
%
%   Applies the coupled Newton Jacobian in real arithmetic. With the
%   holomorphic block A = ML - Jh, the antiholomorphic block B = -Ja, and
%   vc = vr + i*vi, this returns
%       J_real * [vr; vi] = [Re(A*vc + B*conj(vc)); Im(A*vc + B*conj(vc))].
%
%   Inputs:
%     v_real  — stacked [real(vc); imag(vc)]
%     A       — holomorphic Jacobian block ML - Jh
%     B       — antiholomorphic Jacobian block -Ja
%     N       — size of the complex system (length of vc)
%   Outputs:
%     w       — stacked [real(J*vc); imag(J*vc)]
vc = v_real(1:N) + 1i * v_real(N+1:2*N);
wc = A * vc + B * conj(vc);
w  = [real(wc); imag(wc)];
end
