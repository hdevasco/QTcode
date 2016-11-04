function [ minVariance, minPhase, maxVariance, maxPhase ] = optimize_variance( state, initialPhaseMin, initialPhaseMax, S )
%the minimum and maximum quadrature variance of a quantum state
%
%   [ minVariance, minPhase, maxVariance, maxPhase ] = optimize_variance(
%   state, initialPhaseMin, initialPhaseMax )  will find the minimum and
%   maximum quadrature variance (minVariance and max Varaiance)of the
%   state.  It will also return the phase at which the minimum and maximum
%   were found (minPhase and maxPhase).  state may be a column vector or
%   density matrix in Fock basis.  The optimizations will begin with phases
%   initialPhaseMin and initialPhaseMax.  Optional input S is the structure
%   made by init_tables.

if ~exist('S', 'var')
    photons=length(state)-1;
    S = init_tables(photons);
end

options = optimset('Display','off','LargeScale','off','HessUpdate','dfp');

% function that will be minimized
    function mi = minimizeMe(testPhase)
        mi = quadrature_variance(state, testPhase, S);
    end

[minPhase, minVariance] = fminsearch(@minimizeMe, initialPhaseMin, options);

% function that will be maximized
    function ma = maximizeMe(testPhase)
        ma = -quadrature_variance(state, testPhase, S);
    end

% fminsearch has problems with the following.  I don't know why.
[maxPhase, negMaxVariance] = fminunc(@maximizeMe, initialPhaseMax, options);

maxVariance = -negMaxVariance;

minPhase = mod(real(minPhase),2*pi);
maxPhase = mod(real(maxPhase),2*pi);

end