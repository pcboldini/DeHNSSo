function [eigenvalue,index] = EVfilter(EV,V,y,omega,beta,delta99,Mach,mode_select)
% Function which filters out incorrect eigenvalues from an eigenvalue
% spectrum by checking for an exponential fit
%
% The physical disturbance will decay exponentially along "y". Therefore,
% each eigenfunction is compared to:
%       phi = exp( - i sqrt(alpha_r^2 + beta^2) y )
% where alpha is the eigenvalue corresponding to each eigenfunction
% (real component only) and beta the spanwise wavenumber used as input
% for the OS solver.
%
% It is required that the OS solver used a mapping that includes the full
% exponential decay of the eigenfunction! Therefore, make sure that y_max
% of the mapping is sufficiently high.
%
% Optional mode_select:
%   'auto' (default) — minimum-error eigenvalue (no branch preference)
%   'slow' — prefer slow mode S (phase speed c_r = omega/Re(alpha) < 1)
%   'fast' — prefer fast mode F (phase speed c_r > 1)
% Slow/fast filtering: among candidates passing the physical filter, restrict
% to the requested branch by sign of (c_r - 1), then take minimum error.
% Useful at Ma>1 where two discrete modes (Mack F, S) coexist.
%
% M. Kotsonis, K. Groot & J.Y. Boersma - 2017

if nargin<7
    Mach = 0;
end
if nargin<8 || isempty(mode_select)
    mode_select = 'auto';
end

power1 = 6;     % Tweaking factor for the slope error
power2 = 2;     % Tweaking factor for the switch error
treshold = 14;  % logarithm of accuracy of the numerical eigenfunction

use = ~isinf(EV);
EV = EV(use);
V = V(:,use);

NC = length(y);     % length of Chebychev matrix
NE = length(EV);    % length of eigenspectrum

% Parse eigenfunctions V
[Vmax,~] = max(abs(V)); % Store maxima for later use
V = V./repmat(Vmax,NC,1); % Normalize
V(isnan(V))=0+0*sqrt(-1);% replace NaN with zeros

% Combined wavenumbers of alpha and beta
% wavenumber = sqrt(real(EV).^2 + beta^2 - (real(EV)-omega).^2 * Mach)';
wavenumber = real(sqrt(EV.^2 + beta^2 - (EV-omega).^2 * Mach^2))';
%wavenumber = real(sqrt(EV.^2 + beta^2))';
% the former result is derived by Balakumar & Malik (1992) - Discrete modes
% and continuous spectra in supersonic boundary layers, eqs. (3.7) & (2.19)
% note that it is assumed that Pr = 3/4 and Re >> 1, in other cases,
% implement equation (2.22)

% Loop over eigenvalues (and functions)
Error = inf(1,NE);
SlopeError = zeros(1,NE);
Switches = zeros(1,NE);
y_top = zeros(1,NE);
y_bulge = zeros(1,NE);

% ---- Continuous-spectrum branch locations (Balakumar & Malik 1992, eq 3.8) ----
% At finite Ma there are 4 upper-half branch cuts whose branch points are:
%   alpha_vort = omega              (vorticity wave, phase speed = 1)
%   alpha_entr = omega              (entropy wave, degenerate with vorticity at high Re)
%   alpha_acf  = M*omega/(M+1)      (fast acoustic, phase speed = 1 + 1/M)
%   alpha_acs  = M*omega/(M-1)      (slow acoustic, phase speed = 1 - 1/M)
% Discrete modes in a well-resolved LST spectrum can sit close to these; reject
% candidates within a narrow tolerance around each branch.
branch_tol  = 0.02;                 % relative tolerance (2% of omega)
branch_list = omega;                % vorticity/entropy
if Mach > 0
    branch_list = [branch_list, Mach*omega/(Mach+1), Mach*omega/(Mach-1)];
end

for n = 1:NE
    % Skip clearly unphysical eigenvalues early
    if real(EV(n)) < 0 || abs(real(EV(n))) > 1.0 || wavenumber(n) < 1e-10
        continue
    end
    % Reject eigenvalues lying on the continuous-spectrum branches
    if any(abs(real(EV(n)) - branch_list) < branch_tol*abs(omega)) && ...
       abs(imag(EV(n))) < branch_tol*abs(omega)
        continue
    end

    %% Limits in y

    % Determine at which y the numerical accuracy of the eigenfunction will no
    % longer produce an exponential trendline
    y_top(n) = treshold/wavenumber(n);

    % Determine the point where the exponential decay starts and the bulge ends
    % (this point is not strictly defined, just an arbitrary height)
    y_bulge(n) = delta99 + 0.15*delta99/wavenumber(n);
    
    % Check the values
    if y_bulge(n) > y(1)/2
        y_bulge(n) = y(1)/2;
    end
    if y_top(n) < y_bulge(n)*2
        y_top(n) = inf;
    end


    %% Fit of the exponential trend
    % Exponential decay should be equal to the wavenumber

    % Selection "s1" selects the freestream
    s1 = max([find(y<=y_top(n),  1),        find(y<=max(y)*0.85,1)]) ...    % Always skip upper 15%
       : min([find(y>=y_bulge(n),1,'last'), NC-1]);                         % Always skip last node
    EF = log(abs(V(s1,n)));
   
    % Per-segment comparison
    SlopeDifference = wavenumber(n) + diff(EF)./diff(y(s1));
    SegmentLength = -diff(y(s1)); % Weighting factor
    SlopeError_Segment = mean(abs(SegmentLength.*SlopeDifference.^2));

    % Linear least-square fit to determine the slope
    EF_fit = polyfit(y(s1),EF,1);
    SlopeError_LS = abs(wavenumber(n)+EF_fit(1));
    
    % Combine segment comparison and least-square fit
    SlopeError(n) = sqrt(SlopeError_Segment.*SlopeError_LS).^(power1);

    %% Number of switches in gradient

    % Selection "s2" selects the bulge
    s2 = find(y>=y_bulge(n),1,'last'):NC;

    % Count switches, reduced by expected values for each component
    S_imag = countSwitches(imag(V(s2,n)),2);
    S_real = countSwitches(real(V(s2,n)),1);
    S_abs  = countSwitches( abs(V(s2,n)),1);

    % The sum of the number of switches should be zero
    Switches(n) = S_imag+S_real+S_abs;

    % Scale the slope error with the number of switches, using an exponent 
    % because of the logarithmic scaling.
    % Times "power" amplifies the error of multiple switches.
    Error(n) = SlopeError(n).*10.^(Switches(n)*power2);
    

end

%% Select eigenvalue

% Optional slow/fast branch restriction (Ma>1: Mack F-mode vs S-mode).
% Discrete modes are bounded by the continuous-spectrum branches:
%   vorticity:        c_r = 1                    -> Re(alpha) = omega
%   slow acoustic:    c_r = 1 - 1/M              -> Re(alpha) = omega*M/(M-1)
%   fast acoustic:    c_r = 1 + 1/M              -> Re(alpha) = omega*M/(M+1)
% so admissible windows are:
%   slow mode (S):  omega < Re(alpha) < omega*M/(M-1)   [c_r in (1-1/M, 1)]
%   fast mode (F):  omega*M/(M+1) < Re(alpha) < omega   [c_r in (1, 1+1/M)]
% Anything outside is continuous-spectrum or numerical garbage.
switch lower(mode_select)
    case 'slow'
        if Mach > 1
            slow_acf_ub = omega * Mach / (Mach - 1);
            wrong_branch = real(EV) <= omega | real(EV) >= slow_acf_ub;
        else
            wrong_branch = real(EV) <= omega;
        end
        Error(wrong_branch) = inf;
    case 'fast'
        if Mach > 1
            fast_acf_lb = omega * Mach / (Mach + 1);
            wrong_branch = real(EV) >= omega | real(EV) <= fast_acf_lb;
        else
            wrong_branch = real(EV) >= omega;
        end
        Error(wrong_branch) = inf;
    case 'auto'
        % no extra filter
    otherwise
        warning('EVfilter:UnknownModeSelect', ...
                'mode_select=%s not recognised; using ''auto''', mode_select);
end

% Minimum error is the eigenvalue to keep
[~,keep] = min(Error);
if ~isfinite(Error(keep))
    warning('EVfilter:NoCandidate', ...
            'mode_select=%s left no admissible candidate; falling back to all eigenvalues', mode_select);
    % Re-run selection without the branch restriction as a safety net.
    Error_full = Error;  % already-inf'd; rebuild from scratch by recomputing
    % Simplest fallback: take whichever candidate has smallest |Re(EV) - omega|
    [~,keep] = min(abs(real(EV) - omega));
end
eigenvalue = EV(keep);

% Select correct index of original EV list
remain = find(use);
index = remain(keep);


end

function switches = countSwitches(data,expected)
% Count the number of 'switches' in the gradient of 'data', reduced by the
% number of 'expected' switches, up to a minimum of zero.

% Sign of difference of the data segments (hence, removes one element)
sn = sign(diff(data));
% Product of succeeding data segments (also removes one element)
count = sn(1:end-1,:).*sn(2:end,:); % -1 if switched, 1 if not
% Count number of switches
switches = sum(count+1==0);
% Reduce number of switches
switches = abs(switches-expected);
% % Apply minimum of zero
% switches(switches<0) = 0;
end
