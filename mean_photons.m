function [ n ] = mean_photons( rho )
%MEAN_PHOTONS returns the mean number of photons input state
%   Input a state rho, which may be a pure state column vector or density 
%   matrix in Fock
%   basis.  n = mean number of photons in rho.

max_photon = size(rho,1)-1;

if isvector(rho)
    %probabilities will be a column vector of photon number probabilities
    probabilities = rho.*conj(rho);
else
    probabilities = diag(rho);
end

%photon numbers is a row vector [0,1,2,...,max_photon]
photon_numbers = 0:max_photon;
n = photon_numbers*probabilities;
