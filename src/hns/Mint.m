function [Nmat, Mmat, Modevec,Mvec, Nvec] = Mint(M,N)
%MINT  Set up the mode-interaction indexing for harmonic balancing.
%
%   [Nmat, Mmat, Modevec, Mvec, Nvec] = Mint(M, N)
%
%   Builds the tables that record which mode each pairwise interaction feeds.
%   Multiplying two harmonics adds their wavenumbers, so the recipient mode of
%   the interaction between counters i and j carries spectral indices
%   (Mvec(i)+Mvec(j), Nvec(i)+Nvec(j)); Mmat and Nmat store those sums for
%   every pair.
%
%   Inputs:
%     M        — modal truncation in t
%     N        — modal truncation in z
%   Outputs:
%     Nmat     — n-part of the mode affected by each interaction
%     Mmat     — m-part of the mode affected by each interaction
%     Modevec  — [m; n] for every mode counter l
%     Mvec     — m of every mode counter l
%     Nvec     — n of every mode counter l

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
 Nvec((k-1)*(2*M+1)+1:k*(2*M+1)) = -N+(k-1)*ones(1,2*M+1);
end
Modevec = [Mvec;Nvec];

%% Recipient mode of each interaction (wavenumbers add)
for i = 1:L
    for j = 1:L
        Nmat(i,j) = Nvec(i)+Nvec(j);
    end
end
for i = 1:L
    for j = 1:L
        Mmat(i,j) = Mvec(i)+Mvec(j);
    end
end
end
