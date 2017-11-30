clc;
clear all;
% This code estimates the average number of photons from equation
% n_bar = 1/N*(sum(x.^2-1/2));
% where N is the number of quadrature measurements and x is the vector with the quadrature measurements.


% maximum number of photons
 maxPhotonNumber = 10;
 numMeasurements = 20000;
 numAngles       =20;

 angles = pi*(0:numAngles-1)/numAngles;
 angles = repmat(angles,1, ceil(numMeasurements/numAngles));
 angles = angles(1:numMeasurements)';
            
 S     = init_tables(maxPhotonNumber);
 alpha = 1;  % amplitude of coherent states in the superposition
 phase = 0;  % phase between superposition
 psi = generate_cat_vector(alpha, phase, S);

 etaDetector     = 0.9
 etaState        = 0.8;

 % density matrix
 Rho = apply_loss(psi,etaState,S);
 
 Samples = homodyne_samples(-7,7,etaDetector,angles,Rho,S);
 
 % average number photons
 n_bar = n_quadrature(Samples,numMeasurements);