function Measurements = make_measurement_struct(measurementArray, eta, S )
%calculates the POVMs for each measurement and returns a structre
%   make_measurement_struct(measurementArray, eta, S) Input arguments:
%   measurementArray N-by-2, N-by-3, or N-by-4 array whose rows contain
%   (measurement setting, measurement result), (measurement setting,
%   measurement result, number of observations), or (measurement setting,
%   measurement result 1, measurement result 2, number of observations). In
%   the two column case, each measurement is assumed to occur once. In the
%   4 column case, the measurement outcome is described by two numbers (for
%   example the edges of a histogram bin into). eta is the efficiency of
%   the homodyne detector. S = init_tables(max_photons).
%
%   The output structure contains
%   Measurements.nMeasurements = number of different measurement results
%   Measurements.nTotalMeasurements = total number of all measurements
%   Measurements.measurementArray = measurementArray
%   Measurements.counts = list of number of times each outcome is observed
%   Measuremnets.outcomes = array of measurement outcomes
%   Measurements.settings = list of measurment settings
%   Measurements.povmArray = an array of POVMs corresponding to
%                            each measurement result.

nMeasurements = size(measurementArray, 1);
nColumns = size(measurementArray, 2);
settings = measurementArray(:,1);
% If measurementArray has only two columns, assume all counts are "1".
if nColumns==2,
    counts = ones(nMeasurements,1);
    nColumns=3; % tricky.  We do not actually add a column.
else
    counts = measurementArray(:,end);
end
nTotalMeasurements = sum(counts);
outcomes = measurementArray(:,2:(nColumns-1));
nOutCols = size(outcomes,2);

% loss_operation will make the Krauss operators, which we use to correct
% the POVMs for homodyne inefficiency.  These are the Hermitian conjugates
% of the operator elements which one would apply to a state.
if eta == 1
    eHC = 1;
else
    eHC = loss_operation(eta, S, 'ToPOVM');
end

povmArray = zeros(S.dimHilbertSpace,S.dimHilbertSpace,nMeasurements);
if  nOutCols==1,
    for n = 1:nMeasurements
        povmArray(:,:,n) = homodyne_loss_measurement(...
            [settings(n),outcomes(n)], eHC, S, 'return matrix');
    end
elseif nOutCols==2,
    for n = 1:nMeasurements
        povmArray(:,:,n) = coarse_measurement(...
            [settings(n),outcomes(n,1),outcomes(n,2)], eHC, S,...
            'return matrix');
    end
end


Measurements = struct('nMeasurements', nMeasurements, ...
                      'nTotalMeasurements', nTotalMeasurements, ...
                      'settings', settings,...
                      'outcomes', outcomes, ...
                      'counts', counts, ...
                      'povmArray', povmArray);

end