% Calculation of the variance associated with the mean number of photons in the state

alpha =      1;

etaState           =                 0.8;

maxPhotonNumber =                    100;

S         = init_tables(maxPhotonNumber);

n = 0:maxPhotonNumber;

psi = generate_cat_vector(1, 0, 100)

rho = apply_loss(psi,etaState,S);

    
    % p_n constructs the probabilities of p (n).
    p_n = diag(rho);
    
    
 % mean number of photons   
average_numberphotons= mean_photons(rho);

% Variace associated with the mean number of photons
variance =(n-average_numberphotons).^2 *p_n;

% Standard deviation associated with the mean number of photons
sigma = sqrt(variance);

