function d = duplication_antisymmetric(n)
% returns matrix for duplication of a vector into an antisymmetric matrix
%
% This is an adaptation of returns Magnus and Neudecker's duplication 
% matrix of size n for duplicating antisymmetric matrices.  n is the
% dimension of the matrix that is duplicated.
% duplication_antisymmetric(n)*vech(m) = vec(m)

% adapted from code by Thomas P Minka (tpminka@media.mit.edu)

a = tril(ones(n));
i = find(a);
a(i) = 1:length(i);
t = tril(a, -1);
a = a - t';

j = a(:);

m = n*(n+1)/2;
d = zeros(n*n,m);
for r = 1:size(d,1)
    if j(r) > 0
        d(r, j(r)) = 1;
    elseif j(r) < 0
        d(r,abs(j(r))) = -1;
    end
end
