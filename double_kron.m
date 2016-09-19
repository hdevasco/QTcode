function [rho, S] = double_kron( rho1, rho2 )
%double_kron expands Hilbert space and performs kron
%   [rho, S] = double_kron(rho1, rho2)
%   
%   rho = double_kron(rho1, rho2) will double the maximum number of photons
%   in input states rho1 and rho2, and then will perform the kronecker
%   product on them.  rho1 and rho2 may be pure state vectors or density
%   matrices, but they must have the same maximum photon number.  If both
%   rho1 and rho2 are vectors, the output will be a vector, otherwise it
%   will be a density matrix.
%
%   [rho, S] = ... will also output a new S of init_tables for the larger
%   maximum photon number.

photons = length(rho1)-1;

% expand photon number
photons = 2*photons;
dimHilbertSpace = photons+1;
if isvector(rho1)
    rho1(dimHilbertSpace) = 0;
else
    rho1(dimHilbertSpace,dimHilbertSpace) = 0;
end
if isvector(rho2)
    rho2(dimHilbertSpace) = 0;
else
    rho2(dimHilbertSpace,dimHilbertSpace) = 0;
end

% make vectors into density matrices if necessary
if isvector(rho1) && ~isvector(rho2)
    rho1 = rho1*rho1';
end
if isvector(rho2) && ~isvector(rho1)
    rho2 = rho2*rho2';
end

rho = kron(rho1,rho2);
 
if nargout == 2
    S = init_tables(photons);
end

end

