function f = wigner(rho, q, p)
% this function takes a density matrix rho (in the Fock basis) and a set of
% values of quadrature position and momentum q and p (q and p should be the
% same size) and calculates the Wigner function over the q and p using the
% expressions in Leonhardt, _Measuring the Quantum State of Light_, p. 130.

f = zeros(size(q));
r = sqrt( (q.^2) + (p.^2) );
phi = atan2(p, q);
M = size(rho, 1) - 1;

% make tables of the reciprocals of square roots of factorials and square
% roots of n*(n+k) for use by wig()
sqrt_fact_table = zeros(1, M+1);
for k = 0:M
    sqrt_fact_table(k+1) = 1 / prod(sqrt(1:k));
end

n_tab = repmat((0:M).', 1, M+1);
k_tab = repmat(0:M, M+1, 1);
sqrt_table = sqrt(n_tab .* (n_tab + k_tab));

% if input is a pure state vector, convert to density matrix
if isvector(rho)
    rho=rho*rho';
end

for a = 1:size(q,1)
    for b = 1:size(q,2)
        wig_sum = 0;
        for k = (-M):M
            wig_sum = wig_sum + wig(r(a,b), k, rho, M, sqrt_fact_table, sqrt_table) * exp(-1 * 1i * k * phi(a,b));
        end
        f(a,b) = wig_sum;
    end
end

%occasionally, there are small (~10^-16) imaginary parts of f remaining, so we remove
%those

f = real(f);

%TODO: improve speed by removing for loops
%TODO: give faster algorithm for pure states