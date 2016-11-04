function X = binomial_coeffs_matrix(max_n)
% Calculates the binomial coefficients C(i,j) (i,j = 0,...,max_n) using the
% well-known identity C(n,m) = C(n-1, m) + C(n-1, m-1). This is used in the
% calculation of the homodyne measurement operators.

% TODO: try using pascal function to improve speed.

X = zeros(max_n + 1, max_n + 1);
coeffs = zeros(1, max_n + 1);
coeffs(1) = 1;

X(1,:) = coeffs;

for i = 1:max_n
    coeffs = coeffs + [0, coeffs(1:max_n)];
    X(i+1,:) = coeffs;
end