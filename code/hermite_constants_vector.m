function X = hermite_constants_vector(max_n)
% Returns a row vector containing the constant terms of the Hermite
% polynomials for n=0,1,...,max_n.

% it turns out that H_(2n)(0) = (-1)^n * (2n)!/n! and H_(2n+1)(0) = 0...

X = zeros(1, max_n+1);

for n = 0:floor(max_n/2)
    X(2*n + 1) = (-1) ^ n * prod((n+1):(2*n));
end