function w = jacobian_apply(v_real, A, B, N)
%JACOBIAN_APPLY  Newton-Jacobian matrix-vector product in real arithmetic.
%
%   w = jacobian_apply(v_real, A, B, N)
%
%   A = ML - Jh and B = -Ja are the paper's J_h and J_a, built from the nlt_jacobian
%   blocks. With vc = vr + 1i*vi this returns
%       [Re(A*vc + B*conj(vc)); Im(A*vc + B*conj(vc))],
%   which is the real Jacobian [Re(A+B), -Im(A-B); Im(A+B), Re(A-B)] applied to
%   [vr; vi]. v_real stacks [real(vc); imag(vc)]; N is the length of vc.
vc = v_real(1:N) + 1i * v_real(N+1:2*N);
wc = A * vc + B * conj(vc);
w  = [real(wc); imag(wc)];
end
