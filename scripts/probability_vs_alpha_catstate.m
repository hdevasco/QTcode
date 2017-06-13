% Probability of having 10 photons in the Hilbert space as we vary the number of photons in the cat state

P_10 = zeros(9,1);

alpha = [0, 1, 2, 3, 4, 5, 6, 7, 8];

for j= 1:1:9;

psi = generate_cat_vector(j-1, 0, 100)

rho = psi*psi';

w = length(psi);
p_n = zeros(w,1);

for i=1:w-1,
    n = zeros(w,1);
    n(i,1) = 1;
    
    p_n(i,1) = n'*rho*n;
end

% Probability of obtaining 10 photons

P_10(j)= 1-sum(p_n(1:11));

end

plot(alpha,P_10)
