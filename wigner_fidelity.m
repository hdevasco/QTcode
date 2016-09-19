function [ f ] = wigner_fidelity( w1, w2, stepSize )
%wigner_fidelity calculates the overlap integral of two wigner functions
%   wigner_fidelity(w1, w2, stepSize) calculates the overlap of wigner
%   functions w1 and w2.  w1 and w2 are arrays of wigner function values
%   for two quantum states.  The wigner function values arranged on a grid
%   with stepSize distance between neighboring points. Be warned that the
%   overlap integral is equal to the fidelity if w1 or w2 is a pure state.
%   If you know the density matrices for these states, use the function
%   fidelity to compute fidelity.

f = 2*pi*sum(sum(w1.*w2))*stepSize.^2;

end

