function p = purity( state )
%calculates the purity of a quantum state
%   purity(state) calculates the purity of the state state, which is
%   represented as a density matrix.

% p = trace(state^2);

p=trace_of_product(state,state);

end

