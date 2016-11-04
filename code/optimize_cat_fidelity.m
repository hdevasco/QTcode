function [maxFidelity, bestAlpha,bestTheta] = optimize_cat_fidelity(match_state, initial_alpha, initial_phase)
% Finds a pure cat state of type type that has maximum fidelity with input
% density matrix match_state
% input: match_state = density matrix of state you want to match.
%        initial_alpha = initial guess for cat state amplitude
%        If initial_phase = 'even', then the cat is assumed to have phase
%        0, and the phase is not optimized.  If initial_phase = 'odd', then
%        the cat is assumed to have phase pi, and the phase is not
%        optimized.  If initial_phase is a numeric scalar it is used as the
%        initial point for optimizing the phase.
% output: maxFidelity = closest matching fidelity
%         bestAlpha = amplitude of closest cat
%         bestTheta = phase of closest cat
% This function has problems with the minimization, especially when the
% phase is adjustable.

max_photon = length(match_state)-1;
alphaRo = real(initial_alpha);
alphaIo = imag(initial_alpha);

options = optimset('Display','off','LargeScale','off','HessUpdate','dfp');

if ischar(initial_phase)
    if strcmp(initial_phase, 'even')
        theta = 0;
    elseif strcmp(initial_phase, 'odd')
        theta = pi;
    else
        error('Tomography:argChk', 'initial phase string is not even or odd')
    end
    minimizeMe = @(catParameters) -fidelity(match_state,generate_cat_vector(catParameters(1)+catParameters(2)*1i,theta,max_photon));
    [catParameters, minNegFidelity] = fminunc(minimizeMe,[alphaRo,alphaIo],options);
    bestTheta = theta;
elseif isnumeric(initial_phase)
    minimizeMe = @(catParameters) -fidelity(match_state,generate_cat_vector(catParameters(1)+catParameters(2)*1i,catParameters(3),max_photon));
    [catParameters, minNegFidelity] = fminunc(minimizeMe,[alphaRo,alphaIo,initial_phase],options);
    bestTheta = catPararmeters(3);
end

maxFidelity = -minNegFidelity;
bestAlpha = catParameters(1)+1i*catParameters(2);

% TODO: try to replace fminunc with fminsearch so that people without
% optimization toobox can use this.



