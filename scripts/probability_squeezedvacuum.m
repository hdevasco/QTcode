%  Probability of having more than "n" photons in the vacuum state squeezed

% generate state squeezed
psi = generate_squeezed_vacuum_vector(3/4, 100, 'ratio');


rho = psi*psi';

w = length(psi);
p_n = zeros(w,1);

for i=1:w-1,
    n = zeros(w,1);
    n(i,1) = 1;
    
    p_n(i,1) = n'*rho*n;
end

% Probability of the squeezed state to have 10 our more photons 
P_10= 1-sum(p_n(1:11));


% Average number of photons

average_number_photons = mean_photons(rho);


