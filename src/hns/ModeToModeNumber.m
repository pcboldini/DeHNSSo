function [mode] = ModeToModeNumber(m,n,M,N)
%MODETOMODENUMBER  Mode counter for a given (m,n) under truncation (M,N).
%
%   mode = ModeToModeNumber(m, n, M, N)
%
%   Returns the linear mode counter l for the harmonic with spectral indices
%   (m in t, n in z), using the same mode ordering as Mint. The temporal index
%   varies fastest, so l runs over m = -M:M within each spanwise block n.
%
%   Inputs:
%     m     — mode number in t (temporal)
%     n     — mode number in z (spanwise)
%     M     — modal truncation in t
%     N     — modal truncation in z
%   Outputs:
%     mode  — linear mode counter l for (m,n)
%
%   The argument order is (m,n,M,N), matching every call site and the (m,n)
%   labelling used throughout the solver. It previously read (n,m,N,M), which
%   agrees only when M or N is zero and silently returned the wrong counter for
%   runs with both a temporal and a spanwise harmonic.

%% Error check
if abs(m)>M || abs(n)>N
    error('Input exceeds truncation (m,n) > (M,N)')
end

%% Mode counters
L = (2*N+1)*(2*M+1);

% Basic building block in t
M1 = -M:1:M;

% omega (t) mode of each counter l
for i = 1:L
 Mvec(i) = M1(mod(i-1,length(M1))+1);
end

% beta (z) mode of each counter l
for k = 1:L/(2*M+1)
 Nvec([(k-1)*(2*M+1)+1:k*(2*M+1)]) = -N+(k-1)*ones(1,2*M+1);
end

%% Locate the counter matching (m,n)
[~,mindex] = find(Mvec==m);
[~,nindex] = find(Nvec==n);
mode = intersect(mindex,nindex);

end
