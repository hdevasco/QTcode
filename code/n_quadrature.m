function n_bar = n_quadrature(Samples,numMeasurements)
% n_quadrature estimates the average number of photons from the quadrature measurements.
x = Samples(:,2);

n_bar = 1/numMeasurements*(sum(x.^2-1/2));

end