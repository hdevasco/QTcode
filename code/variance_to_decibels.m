function [ decibels ] = variance_to_decibels( variance )
%converts variance of squeezed state to decibels relative to vacuum
%   variance_to_decibels(variance) will convert the variance of a squeezed
%   state to number of decibels of squeezing relative to vacuum.  Vacuum
%   variance is 1/2.  variance should be a real scalar.

vacuumVariance = 1/2;
decibels = 10.*log10(variance./vacuumVariance);

end

