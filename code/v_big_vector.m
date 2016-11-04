function v = v_big_vector(RI)
%vector v used in regularized quadratic approximation of loglikelihood
%   v_big_vector(RI) calculates vector v used in the regularized quadratic
%   approximation of loglikelihood: loglikelihood(rho) =
%   loglikelihood(rho0)+v.'*aRI + second order terms.  RI is a structre
%   containing data derived from a density matrix rho0 and R matrix.

v = 2*(-RI.traceRRho0*RI.rhoSqrtRI+RI.qRI);
v = real(v);

end