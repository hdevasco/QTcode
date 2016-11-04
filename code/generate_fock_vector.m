function psi = generate_fock_vector(n, maxPhoton)
% Creates photon number state
%   psi = generate_fock_vector(n, maxPhoton) returns the pure state vector
%   for photon number eigenstate n in the photon number basis. maxPhoton is
%   the photon number at which the Hilbert space is truncated.

psi = zeros(maxPhoton+1,1);
psi(n+1,1) = 1;