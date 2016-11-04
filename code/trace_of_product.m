function [ tp ] = trace_of_product( A, B )
%TRACE_OF_PRODUCT calculates that trace of the product of two matrices.
%   trace_of_product calculates trace(A*B), but it is faster than
%   trace(A*B) for large matrices.  I'm not sure where the cut off is, but
%   for small matrices trace(A*B) seems to be faster.
%   Input: A and B are 2D arrays which can be multiplied.
%   Output: Trace(A*B)

tp = sum(sum(A.*B.'));
    
end

