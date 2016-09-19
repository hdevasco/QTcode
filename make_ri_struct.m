function RI = make_ri_struct( R, Measurements, makeNoM )
%makes a structure using current density matrix, for rga_iterations
%   make_ri_struct( R, Measurements) makes a structure of data derived from density
%   matrix rho and matrix R.  R should be the structure made by
%   make_r_struct.  Measurements is the structure of measurement data.
%   RI is used for several functions in the regularized
%   gradient ascent family.

rhoRI = vectorize_r_i(R.rho);
traceRRho0 = Measurements.nTotalMeasurements; % previously was traceRRho0 = trace_of_product(R.r, R.rho);
rhoSqrt = sqrtm(R.rho);
rhoSqrtR = real(rhoSqrt);
rhoSqrtI = imag(rhoSqrt);
rhoSqrtRI = vectorize_r_i(rhoSqrt);
qRI = vectorize_r_i(R.r*rhoSqrt);

RI = struct('rhoRI', rhoRI, ...
            'traceRRho0', traceRRho0, ...
            'rhoSqrt', rhoSqrt, ...
            'rhoSqrtR', rhoSqrtR, ...
            'rhoSqrtI', rhoSqrtI, ...
            'rhoSqrtRI', rhoSqrtRI, ...
            'qRI', qRI);

v = v_big_vector(RI);

if exist('makeNoM','var')
    m=zeros(max(size(v)));    
else
    m = m_big_matrix(R, RI, Measurements);
end   
 
[mDiagonalizer, mDiag] = eig(m);
mEigList = diag(mDiag);
mMaxEigenvalue = max(mEigList);

if mMaxEigenvalue == 0
    mMaxEigenvalue = 1;
end

RI.v = v;
RI.m = m;
RI.mDiagonalizer = mDiagonalizer;
RI.mEigList = mEigList;
RI.mMaxEigenvalue = mMaxEigenvalue;

end

