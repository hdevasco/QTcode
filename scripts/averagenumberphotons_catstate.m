% (<Psi|n|Psi>) Average number of protons according to the alpha value (catstate)

average_number_photons = zeros(9,1);

alpha = [0,1, 2, 3, 4, 5, 6, 7, 8];

for j=1:1:9;

psi = generate_cat_vector(j-1, 0, 100)

rho = psi*psi';

average_number_photons(j)= mean_photons(rho);

end

plot(alpha, average_number_photons)