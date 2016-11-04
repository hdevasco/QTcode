function maxEig = max_eig( matrix )
% computes the maximum eigenvalue of a matrix
%   max_eig(matrix) computes the maximum value of the matrix.

eigenvalues = eig(matrix);
maxEig = max(eigenvalues);

end

