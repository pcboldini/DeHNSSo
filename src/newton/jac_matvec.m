function w = jac_matvec(v_real, A, B, N)
%JAC_MATVEC  Real-arithmetic matvec for the HNS Newton Jacobian.
%
%   w = jac_matvec(v_real, A, B, N)
%
%   Applies the coupled Newton Jacobian in real arithmetic. With the
%   holomorphic part A = ML - Jh and antiholomorphic part B = -Ja, and
%   vc = vr + i*vi, returns
%       J_real * [vr; vi] = [Re(A*vc + B*conj(vc)); Im(A*vc + B*conj(vc))].
%
%   Inputs:
%     v_real  — stacked [real(vc); imag(vc)]
%     A       — holomorphic Jacobian block ML - Jh
%     B       — antiholomorphic Jacobian block -Ja
%     N       — complex system size (length of vc)
%   Outputs:
%     w       — stacked [real(J*vc); imag(J*vc)]
vc = v_real(1:N) + 1i * v_real(N+1:2*N);
wc = A * vc + B * conj(vc);
w  = [real(wc); imag(wc)];
end
