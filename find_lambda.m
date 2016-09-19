function lambda = find_lambda(RI, stepSize);
% finds labmda to provide desired stepSize for RGA
%   lambda = find_lambda(RI, stepSize) finds a lambda to match the
%   desired step size in regularized gradient ascent iterations.
%   It uses Haley's method to accomplish this.  RI is the structure
%   made by make_ri_struct.  stepSize is the desired step size.

lambda = 1.001/2*max([RI.mMaxEigenvalue; 0]);

for n = 1:3
  lm = 2*lambda-RI.mEigList;
  lmInv = 1./lm;
  lmInv2 = lmInv.^2;
  lmInv3 = lmInv.^3;
  lmInv4 = lmInv.^4;
  mDv = RI.mDiagonalizer.'*RI.v;
  mDvT = mDv.'
  curentStep = mDvT*diag(lmInv2)*mDv-stepSize;
  firstDerivative = -2*mDvT*diag(lmInv3)*mDv;
  secondDerivative = 6*mDvT*diag(lmInv4)*mDv;
  lambda = (lambda-2*currentStep*firstDerivative)/...
           (2*firstDerivative^2-currentStep*secondDerivative);
end
