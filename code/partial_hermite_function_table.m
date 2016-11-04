function X = partial_hermite_function_table(max_n, bin_cof, hermite_consts)
% Using a table of binomial coefficients (bin_cof), a list of constant
% terms for Hermite polynomials (hermite_consts), and a maximum n (max_n),
% this function computes a table that can be used to speed up calculating
% the Hermite function
%
% phi_n(x) = exp(-x^2 / 2) * H_n(x) / sqrt(2^n * n! * sqrt(pi)).

X = bin_cof;

for n = 0:max_n
    % 
    X(n+1,:) = X(n+1,:) ./ prod(sqrt(1:n));
    
    X(n+1,1:(n+1)) = X(n+1,1:(n+1)) .* hermite_consts(1:(n+1)) .* 2.^((n/2) - (0:n));
end

X = X .* pi^(-1/4);