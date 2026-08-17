function [ HB ] = Hbalancing(M, N, Mmat, Nmat, Modevec )
%HBALANCING  Build the harmonic-balancing registry for mode interactions.
%
%   HB = Hbalancing(M, N, Mmat, Nmat, Modevec)
%
%   Returns a stack of binary matrices used to assemble the mode-to-mode
%   nonlinear forcing of the HNS solver. In slice HB(:,:,j) a 1 marks a pair
%   of modes (row, column) whose interaction forces mode j. For N = 2, M = 0
%   the slice for the mean-flow distortion (0,0) reads:
%
%       n,m  -2,0  -1,0   0,0   1,0   2,0
%      -2,0    0     0     0     0     1
%      -1,0    0     0     0     1     0
%       0,0    0     0     1     0     0
%       1,0    0     1     0     0     0
%       2,0    1     0     0     0     0
%
%   since (-2,0)x(2,0), (-1,0)x(1,0), ... all force (0,0). In practice the
%   (0,0)x(0,0) self-interaction is excluded from the forcing calculation.
%
%   Inputs:
%     M        — number of modes in t
%     N        — number of modes in z
%     Mmat     — m-part of the modes affected by each interaction
%     Nmat     — n-part of the modes affected by each interaction
%     Modevec  — (m,n) of each mode counter l
%   Outputs:
%     HB       — binary interaction registries stacked over recipient modes

%% Build interaction registry
L = (2*N+1)*(2*M+1);
for i = 1:L
% Modes (m,n) received by mode i (counter l)
m=Modevec(1,i);
n=Modevec(2,i);
A = Mmat == m;      % entries with the matching m
B = Nmat == n;      % entries with the matching n
HB(:,:,i) = A.*B;
end

end
