%  Probability of having more than "n" photons in the vacuum state squeezed

% generate state squeezed
psi = generate_squeezed_vacuum_vector(3/4, 100, 'ratio');

maxPhotonNumber = 100;
S               = init_tables(maxPhotonNumber);

% Apply_loss returns the density matrix after passing through a lossy medium
rho = apply_loss(psi, 0.8, S);

w = length(rho);
p_n = zeros(w,1);


for i=1:w-1,
    % p_n constructs an array whose lines are the probabilities of p (n).
    p_n(:,1) = diag(rho);
    
end

% Probability of the squeezed state to have 10 our more photons 
P_10= 1-sum(p_n(1:11));


% Average number of photons

average_number_photons = mean_photons(rho);


