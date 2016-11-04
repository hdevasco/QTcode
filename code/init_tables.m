function S = init_tables(photons)
% Calculates a number of useful tables and puts them into a structure
% for convenience.

S.photons = photons;
S.dimHilbertSpace = photons+1;
S.binom = binomial_coeffs_matrix(photons);
S.hermite = hermite_constants_vector(photons);
S.part_hf = partial_hermite_function_table(photons, S.binom, S.hermite);
S.h_exp = h_exponents_table(photons);


if photons > 1
% meanQuadratureMatrix is used to calculate the expectation value of
% quadrature x from the photon number density matrix
% <x>=Tr(rho*meanQuadratureMatrix)
expXMatrixElements = realsqrt((1:photons)/2);
S.expXMatrix = diag(expXMatrixElements,1)+diag(expXMatrixElements,-1);

% expXSquaredMatrix is used to calculate the expectation value of x^2
% <x^2>=Tr(rho*expXSquaredMatrix)
expXSquaredDiagonals = (0:photons)+0.5;
expXSquaredOffDiagonals = realsqrt(S.binom(3:(photons+1),3)*0.5);
S.expXSquaredMatrix = diag(expXSquaredDiagonals)+diag(expXSquaredOffDiagonals,2)+diag(expXSquaredOffDiagonals,-2);
end
