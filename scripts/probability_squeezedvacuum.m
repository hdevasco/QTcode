
%  Probability of having more than "n" photons in the vacuum state squeezed

% generate state squeezed
psi = generate_squeezed_vacuum_vector(3/4, 100, 'ratio');

maxPhotonNumber = 100;
S               = init_tables(maxPhotonNumber);

rho = apply_loss(psi, 1/2, S);

w = length(psi);
p_n = zeros(w,1);


for i=1:w-1,
    n = zeros(w,1);
    n(i,1) = 1;
    p_n = abs(psi).^2;

end

%Probability of photons
p_photons = diag(rho);

% Probability of the squeezed state to have 10 our more photons 
P_10= 1-sum(p_n(1:11));


% Average number of photons

average_number_photons = mean_photons(rho);


