function [ boundList ] = confidence_interval( Measurements, photons, eta, confidenceLevel, linearForm, maxIterations, stop, DiagnosticsMl)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if isnumeric(photons)
    S = init_tables(photons);
elseif isstruct(photons)
    S = photons;
end

if isnumeric(Measurements)
    Measurements = make_measurement_struct(Measurements, eta, S);
end

if ~exist('DiagnosticsMl', 'var')
    [rhoMl, DiagnosticsMl] = combined_optimization(Measurements, S, eta, 0, maxIterations, stop);
else
    rhoMl = DiagnosticsMl.rhoArray(:,:,end);
end
maxLogLikelihood = DiagnosticsMl.loglikelihoodList(end);
maxLikelihoodLinearTerm = real(trace_of_product(rhoMl, linearForm));

rhoAlpha = rhoMl;

for alpha = [-1,1]

    alphaLowerBound = 0;
    alphaUpperBound = inf;

    while alphaUpperBound == inf
        [rhoAlpha, DiagnosticsAlpha] = combined_optimization(Measurements, photons, ...
                                        eta, alpha*linearForm, maxIterations, stop, rhoAlpha);
        logLikelihoodAlpha = DiagnosticsAlpha.loglikelihoodList(end)
        pAlpha = 1 - exp(logLikelihoodAlpha-maxLogLikelihood)
        if pAlpha < confidenceLevel
            alphaLowerBound = alpha;
            alpha = 2*alpha;
        elseif pAlpha >= confidenceLevel;
            alphaUpperBound = alpha;
        end
    end

    disp('found bounds')


    while abs(pAlpha-confidenceLevel) > 0.01
        %pAlpha-confidenceLevel
        alpha = 0.5*(alphaUpperBound + alphaLowerBound);
        [rhoAlpha, DiagnosticsAlpha] = combined_optimization(Measurements, ...
                                                          photons, eta, alpha*linearForm, maxIterations, stop, rhoAlpha);
        %alpha
        if DiagnosticsAlpha.linearTermList(end) == 0 && alpha == 22
          save('check_linear_term.mat', 'DiagnosticsAlpha');
        end
        DiagnosticsAlpha.linearTermList(end)/alpha;
        logLikelihoodAlpha = DiagnosticsAlpha.loglikelihoodList(end);
        pAlpha = 1-exp(logLikelihoodAlpha-maxLogLikelihood);
        if pAlpha < confidenceLevel - 0.01
            alphaLowerBound = alpha;
        elseif pAlpha > confidenceLevel + 0.01
            alphaUpperBound = alpha;
        end
    end

    %pAlpha - confidenceLevel;

    if alpha > 0
        upperBound = DiagnosticsAlpha.linearTermList(end)/alpha;
    else
        lowerBound = DiagnosticsAlpha.linearTermList(end)/alpha;
    end
end
boundList = [lowerBound; maxLikelihoodLinearTerm; upperBound];
    
end

