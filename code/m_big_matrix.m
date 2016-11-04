    function m = m_big_matrix(R, RI, Measurements)
%calculates the big matrix m used for the regularized gradient ascent
%   m_big_matrix(R, RI, Measurements ) calculates the big matrix m used for
%   the regularized gradient ascent.  It appears in the quadratic expansion
%   of the loglikelihood function loglik_regularized.  It is the hessian
%   (matrix of second derivatives) with respect to the real and imaginary
%   parameters of the vector a.  R is the structure containing data about
%   the current density matrix and likelihood. RI is a structre containing
%   data derived from density matrix rho0 and matrix R. Measurements is the
%   structre of measurement results and povms.

dimHilbertSpace = size(R.r, 1);
rR = real(R.r);
rI = imag(R.r);

% second order relations between density matrix and aRI.
m1 = 2*(kron(kron(eye(2),eye(dimHilbertSpace)),rR) ...
        +kron(kron([0,-1;1,0],eye(dimHilbertSpace)),rI));
dimASpace = size(RI.rhoSqrtRI,1);
m2 = -2*RI.traceRRho0*eye(dimASpace);
m3 = -8*RI.rhoSqrtRI*RI.qRI.';
m4 = 8*RI.traceRRho0*(RI.rhoSqrtRI*RI.rhoSqrtRI.');
alpha = -Measurements.measurementArray(:,3)./R.tprl.^2;

% calculating the hessian in the space of deltaRI.
% the povm's are hermitian, so I break them into symmetric real and
% antisymmetric complex parts.  I then vectorize the parts and remove
% redundant entries.  The matrices ds and da restore the redundant parts.
t = transposing_vectorized(dimHilbertSpace);
lowerDiagonals = logical(tril(ones(dimHilbertSpace)));
ds = duplication_symmetric(dimHilbertSpace);
da = duplication_antisymmetric(dimHilbertSpace);
hDim = dimHilbertSpace*(dimHilbertSpace+1)/2;
hRr = zeros(hDim);
hRi = zeros(hDim);
hI = zeros(hDim);
for iPovm = 1:Measurements.nMeasurements
    povm = Measurements.povmArray(:,:,iPovm);
    povm = povm(lowerDiagonals); % makes vech(povm)
    povmR = real(povm);
    povmI = imag(povm);
    hRr = hRr + (alpha(iPovm)*povmR)*povmR.';
    hRi = hRi + (alpha(iPovm)*povmI)*povmI.';
    hI = hI + (alpha(iPovm)*povmI)*povmR.';
end
hR = t*(ds*hRr*ds.'-da*hRi*da.')*t;
hI = da*hI*ds.';
hI = t*(hI+hI.')*t;

%  the following prepares the matrix deltaAmatrix.  rhoNew = rhoOld + delta
%  deltaRI = deltaAmatrix*aRI.  In the space parameterized with aRI, the
%  hessian becomes deltaAmatrix.'*hessian*deltaAmatrix.
deltaAmatrix1 = 2*RI.rhoRI*RI.rhoSqrtRI.'*[-t,zeros(dimHilbertSpace^2);t,zeros(dimHilbertSpace^2)];
deltaAmatrix2DiagonalPart = kron(eye(dimHilbertSpace),RI.rhoSqrtR)*t;
deltaAmatrix2OffDiagonalPart = kron(eye(dimHilbertSpace),RI.rhoSqrtI)*t;
deltaAmatrix2 = [deltaAmatrix2DiagonalPart, deltaAmatrix2OffDiagonalPart; deltaAmatrix2OffDiagonalPart, -deltaAmatrix2DiagonalPart];
deltaAmatrix3DiagonalPart = kron(RI.rhoSqrtR.', eye(dimHilbertSpace));
deltaAmatrix3OffDiagonalPart = kron(RI.rhoSqrtI.', eye(dimHilbertSpace));
deltaAmatrix3 = [deltaAmatrix3DiagonalPart,-deltaAmatrix3OffDiagonalPart; deltaAmatrix3OffDiagonalPart, deltaAmatrix3DiagonalPart];
deltaAmatrix = deltaAmatrix1 + deltaAmatrix2 + deltaAmatrix3;

hessian = deltaAmatrix.'*[hR,-hI;-hI,-hR]*deltaAmatrix;






m = m1+m2+m3+m4+hessian;

% m should be symmetric and real
m = real(m);
m = 0.5*(m+m.');

end

% Here is an earlier method to calculate the hessian matrix.
% 
% m5 = zeros(dimASpace);
% for iPovm = 1:Measurements.nMeasurements
%     iJ = Measurements.povmArray(:,:,iPovm)*RI.rhoSqrt;
%     % iJRI = vectorize_r_i(iJ), but writing the code here allows Matlab to
%     % process the loop faster.
%     iJvector = iJ(:);
%     iJR = real(iJvector);
%     iJI = imag(iJvector);
%     iJRI = [iJR;iJI];
%     iVector = 2*(-R.tprl(iPovm)*RI.rhoSqrtRI+iJRI);
%     iTerm = (alpha(iPovm)*iVector)*iVector.';
%     m5 = m5 + iTerm;
% end
% 
% Many awesome tricks to vectorize for loop above!!
% They are actually slower than the for loop, and use far too much memory.
%
%
% hRrV = zeros(hDim);
% hRiV = zeros(hDim);
% hIV = zeros(hDim);
% lowerDiagonals = repmat(lowerDiagonals, [1, 1, Measurements.nMeasurements]);
% povmArray = Measurements.povmArray(lowerDiagonals);
% povmArray = reshape(povmArray, [hDim, 1, Measurements.nMeasurements]);
% povmArrayR = real(povmArray);
% povmArrayRTp = permute(povmArrayR, [2, 1, 3]);
% povmArrayI = imag(povmArray);
% povmArrayITp = permute(povmArrayI, [2, 1, 3]);
% alphaReshape = reshape(alpha, [1,1,Measurements.nMeasurements]);
% povmArrayRAlpha = bsxfun(@times,alphaReshape, povmArrayR);
% povmArrayIAlpha = bsxfun(@times,alphaReshape, povmArrayI);
% hRrV = bsxfun(@times, povmArrayRAlpha, povmArrayRTp);
% hRrV = sum(hRrV, 3);
% hRiV = bsxfun(@times, povmArrayIAlpha, povmArrayITp);
% hRiV = sum(hRiV, 3);
% hIV = bsxfun(@times, povmArrayIAlpha, povmArrayRTp);
% hIV = sum(hIV, 3);
% hRV = t*(ds*hRrV*ds.'-da*hRiV*da.')*t;
% hIV = da*hIV*ds.';
% hIV = t*(hIV+hIV.')*t;