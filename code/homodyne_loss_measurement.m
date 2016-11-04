function x = homodyne_loss_measurement(paramList, eta, S, returnMatrix)
% Calculates the measurement operator for a homodyne measurement
%   x = homodyne_loss_measurement(paramList, eta, S, returnMatrix)
%   calculates the measurement operator for a homodyne measurement. 
%   paramList(1) is the local oscillator's phase and paramList(2) is the
%   measurement result corresponding to measurement result of paramList(2).
%   S is the structure created by init_tables.  eta = efficiency of the 
%   homodyne measurement, or eta can be the 3D array containing the 
%   operation elements used in the operator sum that transforms the
%   measurement through transmissivity eta.  In the later case, eta should
%   be the hermitian conjugate of the operator elements that would
%   transform a state.  If the measurement is a pure state projection, then
%   x will be a vector, unless the option 'return matrix' is used.  Then it
%   returns the matrix measurement operator.  If the measurement is mixed,
%   then it always returns a matrix.

theta = paramList(1);
x = paramList(2);

% calculate the Hermite function H_n(x) for n=0,...,S.maxPhotons
A = S.part_hf .* x .^ S.h_exp;
x = sum(A, 2) .* exp((-1/2) .* (x^2));
x = x .* exp(1i .* theta .* (0:(S.photons)).');
% now x(n) = e^{i*n*theta}*psi_n(x)

if ~isequal(eta,1)
    x = apply_loss(x, eta, S, 'DoHermitianConjugate');
end

if strcmp(returnMatrix,'return matrix')&&isvector(x)
    x = x*x';
end