function [mode] = ModeToModeNumber(n,m,N,M)
%MODETOMODENUMBER  Mode counter for a given (n,m) under truncation (N,M).
%
%   mode = ModeToModeNumber(n, m, N, M)
%
%   Returns the linear mode counter l for the harmonic with spectral indices
%   (m in t, n in z), using the same mode ordering as Mint.
%
%   Inputs:
%     n     — mode number in z
%     m     — mode number in t
%     N     — modal truncation in z
%     M     — modal truncation in t
%   Outputs:
%     mode  — linear mode counter l for (m,n)

%% Error check
if n>N || m>M
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
