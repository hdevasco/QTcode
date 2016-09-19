function [ t ] = transposing_vectorized( d )
%calculates the matrix needed for transposing a vectorized matrix
%   transposing_vectorized(d) generates a matrix with the function
%   t*vec(m)= vec(m.'), where m is a square matrix with dimension d.

dSquared = d^2;
t = zeros(dSquared,dSquared);
for a = 1:dSquared
    for b = 1:dSquared
        if b == 1 + d*(a-1)-(d^2-1)*floor((a-1)/d)
            t(a,b)=1;
        end
    end
end


end

