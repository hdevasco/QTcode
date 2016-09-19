function [ vri ] = vectorize_r_i( matrix )
%vectorize_r_i vectorizes a matrix and separates real and imaginary parts
%   vectorize_r_i(matrix) returns a 2*n*m element vector that contains the
%   real parts linearly indexed matrix followed by the imaginary parts.

matrixVectorized = matrix(:);

vr = real(matrixVectorized);
vi = imag(matrixVectorized);

vri = [vr;vi];

end

