function X = h_exponents_table(max_n)
% Calculates a table of exponents that comes in useful when calculating the
% Hermite function.

X = zeros(max_n+1);

for n = 0:max_n
    X(n+1,1:(n+1)) = n - (0:n);
end