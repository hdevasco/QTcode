psi = generate_squeezed_vacuum_vector(3/4, 10, 'ratio');

rho = psi*psi';

% n = <n>: average number of photons

n = mean_photons(rho);