function sigma = deviation_n_rho(rho)
%Deviation_n_rho returns the standard deviation associated with the mean number 
%of photons of the rho state, and uses the rho matrix as input.

maxPhotonNumber = size(rho,1)-1;

if isvector(rho)
    %probabilities will be a column vector of photon number probabilities
    p_n = rho.*conj(rho);
else
    p_n = diag(rho);
end


 n   = 0:maxPhotonNumber;
 
 %  mean number of photons   
 average_numberphotons= mean_photons(rho);
 
% Variace associated with the mean number of photons
variance =(n-average_numberphotons).^2 *p_n;
% 
% Standard deviation associated with the mean number of photons
 sigma = sqrt(variance);

