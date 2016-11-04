function [ t ] = trace_distance( rho, sigma )
%trace_distance calculates the trace distance between two density matrices
%   trace_distance(rho, sigma) is the trace distance between density matrices
%   rho and sigma.

t = 1/2*trace(abs(rho-sigma));

end

