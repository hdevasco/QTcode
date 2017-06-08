% Probabilidade de termos mais de "n" fótons no estado de vácuo comprimido

%  Probability of having more than "n" photons in the vacuum state squeezed

psi = generate_squeezed_vacuum_vector(3/4, 10, 'ratio');

% generate state squeezed
rho = psi*psi';

w = length(psi);
p_n = zeros(w,1);

for i=1:w-1,
    n = zeros(w,1);
    n(i,1) = 1;
    
    p_n(i,1) = n'*rho*n;
end


P_10= 1-sum(p_n(1:11));