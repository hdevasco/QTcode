function sigma = deviation_n(rho)
%Deviation_n returns the standard deviation associated with the mean number 
%of photons of the rho state, and uses the rho matrix as input.

n   = 0:maxPhotonNumber;

    % p_n constructs the probabilities of p (n).
    p_n = diag(rho);
    
    
 % mean number of photons   
average_numberphotons= mean_photons(rho);

% Variace associated with the mean number of photons
variance =(n-average_numberphotons).^2 *p_n;

% Standard deviation associated with the mean number of photons
sigma = sqrt(variance);
