function Measurements = make_measurement_struct(measurementArray, eta, S )
%calculates the POVMs for each measurement and returns a structre
%   make_measurement_struct(measurementArray, eta, S)
%   Input arguments: measurementArray N-by-2 or N-by-3 array whose rows 
%   contain (measurement setting, measurement result, ...
%   number of observations), but if the third column is not included each
%   measurement is assumed to occur once. eta is the efficiency of the 
%   homodyne detector. S = init_tables(max_photons).
%
%   Measurements.nMeasurements = number of different measurement results
%   Measurements.nTotalMeasurements = total number of all measurements
%   Measurements.measurementArray = measurementArray
%   Measurements.povmArray = an array of POVMs corresponding to
%                            each measurement result.

nMeasurements = size(measurementArray, 1);

% if measurementArray has only two columns, add column of 1s to 
% measurementArray to show each observation occurred once
if size(measurementArray,2)==2;
    measurementArray = [measurementArray,ones(nMeasurements,1)];
end

nTotalMeasurements = sum(measurementArray(:,3));

povmArray = zeros(S.dimHilbertSpace,S.dimHilbertSpace,nMeasurements);

for n = 1:nMeasurements
    povmArray(:,:,n) = homodyne_loss_measurement(measurementArray(n,1:2), ...
                                                 eta, S, 'return matrix');
end

Measurements = struct('nMeasurements', nMeasurements, ...
                      'nTotalMeasurements', nTotalMeasurements, ...
                      'measurementArray', measurementArray, ...
                      'povmArray', povmArray);

end
