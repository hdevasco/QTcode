function aRI = ari_of_lambda(lambda, RI)
% calculates vector aRI that maximizes loglikelihood function
%   ari_of_lambda(lambda, RI) calculates the vector aRI that maximizes the
%   loglikelihood function, subject to the constraint that
%   aRI.'*aRI=stepSize. lambda is the Lagrange multiplier appearing
%   when optimizing loglikelihood function over aRI, subject to the
%   constrat above.  RI is the structure made by make_ri_struct; it
%   contains data related to the current density matrix position expressed
%   in the real parameterized space.

lm = lambda-RI.mEigList;
lmInv = 1./lm;
lmInv = diag(lmInv);

%   v is the vector of first order coefficients,
%   calculated by v_big_vector, and mEigList is the list of eigenvalues of
%   m, and mDiagonalizer is the unitary matrix that diagonalizes m.
%   L(aRI) = L(rho0)+v.'*aRI+1/2*aRI.'*m*aRI.

aRI = RI.mDiagonalizer*lmInv*RI.mDiagonalizer.'*RI.v;