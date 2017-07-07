% Probability of having 10 photons in the Hilbert space as we vary the number of photons in the cat state

P_10 =                        zeros(9,1);
average_numberphotons =       zeros(9,1);

alpha =      [0, 1, 2, 3, 4, 5, 6, 7, 8];

etaState           =                 0.8;

maxPhotonNumber =                    100;
S         = init_tables(maxPhotonNumber);

for j= 1:1:9;

psi = generate_cat_vector(j-1, 0, 100)

rho = apply_loss(psi,etaState,S);

w = length(rho);
p_n = zeros(w,1);

for i=1:w-1,
    
    % p_n constructs an array whose lines are the probabilities of p (n).
    p_n(:,1) = diag(rho);
    
end


% Probability of obtaining 10 our more photons

P_10(j)= 1-sum(p_n(1:11));


average_numberphotons(j)= mean_photons(rho);

end

plot(alpha,P_10)
