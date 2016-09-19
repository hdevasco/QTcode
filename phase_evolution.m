function [ u ] = phase_evolution( theta, photons )
%phase_evolution makes a phase evolution matrix with phase theta
%   phase_evoultion(theta, photons) makes a phase evolution unitary matrix
%   for phase theta.  The matrix will operate in the photons+1 dimensional
%   space.

u = exp(-1i*theta*(0:photons));
u = diag(u);

end

