function w = wigner_origin_matrix(dimHilbertSpace)
%makes matrix to calculate wigner function at origin
%    wigner_origin_matrix(dimHilbertSpace) creates the matrix w,
%    which should be used to calculate the Wigner function at the
%    origin.  W(0,0) = trace(rho*w).

w = [1,-1];
dHalf = dimHilbertSpace/2;
dHalfInt = floor(dHalf);

wDiag = repmat(w,[1,dHalfInt]);
if 2*dHalfInt ~= dimHilbertSpace
  wDiag = [wDiag,1];
end
wDiag = wDiag/pi;
w = diag(wDiag);
