function [ maxFidelity, bestVarianceRatio, bestPhase ] = optimize_squeezed_fidelity( matchState, initialVarianceRatio, initialPhase )
%optimize_squeezed_fidelity finds squeezed state with maximum fidelity
%   optimize_squeezed_fidelity(matchState, initialVarianceRatio, initialPhase)
%   will find the squeezed vacuum state with maximum fidelity with
%   matchState.  matchState should be a density matrix.  initialVarianceRatio and
%   initialPhase are used to initialize the maximization routine.
%   initialVarianceRatio is the level of squeezing, and initialPhase is the phase
%   evolution of the squeezed state.  Output maxFidelity is the fidelity of
%   match_state and the optimial squeezed state, which has squeezing
%   parameter bestZeta and phase bestPhase.

maxPhoton = length(matchState)-1;
maxPhoton =4*maxPhoton;
if isvector(matchState)
    matchState(maxPhoton+1)=0;
else
    matchState(maxPhoton+1,maxPhoton+1) = 0;
end

options = optimset('Display','off','LargeScale','off','HessUpdate','dfp');

% function that will be minimized: negative of the fidelity of matchState
% with a squeezed state
    function m = minimizeMe(testSqueezedParameters)
        testVarianceRatio = testSqueezedParameters(1);
        testPhase = testSqueezedParameters(2);
        squeezedVector = generate_squeezed_vacuum_vector(testVarianceRatio, maxPhoton);
        squeezedVector = phase_evolution(testPhase, maxPhoton)*squeezedVector;
        m = -fidelity(matchState, squeezedVector);
    end

warning off fockSpace:tooSmall

initialSqueezedParameters = [initialVarianceRatio, initialPhase];
[bestSqueezedParameters, minNegFidelity] = fminsearch(@minimizeMe, initialSqueezedParameters, options);

warning on fockSpace:tooSmall

bestVarianceRatio = bestSqueezedParameters(1);
bestPhase = bestSqueezedParameters(2);
maxFidelity = -minNegFidelity;

end