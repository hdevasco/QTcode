function a = unvectorize_r_i( aRI )
%
%   unvectorize_r_i(aRI) takes as input a column vector and returns the
%   complex matrix made of aRI's elements.  This function undoes the work
%   of vectorize_r_i, which breaks a matrix into real and complex parts and
%   vectorizes it.

dimASpace = size(aRI,1);
dimHilbertSpace = sqrt(0.5*dimASpace);

aR = aRI(1:(dimASpace/2));
aI = aRI((dimASpace/2+1):end);
aComplex = aR + 1i*aI;
a = reshape(aComplex, [dimHilbertSpace, dimHilbertSpace]);

end

