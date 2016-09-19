function [ variance ] = decibels_to_variance( decibels )
%converts decibels of squeezing to variance of squeezed state
%   decibels_to_variance(dB) returns variance of squeezed state that has
%   been squeezed to dB decibels relative to the vacuum.  Vacuum variance
%   is 1/2.  dB should be a real scalar.

vacuumVariance = 0.5;

variance = vacuumVariance.*10.^(decibels./10);


end

