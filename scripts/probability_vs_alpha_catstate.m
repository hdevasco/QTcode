% Probability of having 15 photons in the Hilbert space as we vary the number of photons in the cat state

P_15 =                        zeros(9,1);
average_numberphotons =       zeros(9,1);

alpha =      [0, 1, 2, 3, 4, 5, 6, 7, 8];

etaState           =                 0.8;

maxPhotonNumber =                    100;
S         = init_tables(maxPhotonNumber);

for j= 1:1:9;

psi = generate_cat_vector(j-1, 0, 100)

rho = apply_loss(psi,etaState,S);


    
    % p_n constructs the probabilities of p (n).
    p_n = diag(rho);
    


% Probability of obtaining 10 our more photons

% P_10(j)= 1-sum(p_n(1:11));

P_15(j)= 1-sum(p_n(1:16));


average_numberphotons(j)= mean_photons(rho);

end

plot(alpha,P_15)
