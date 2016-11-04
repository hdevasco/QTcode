function v = quadrature_variance( state, phase, S )
%the quadrature variance of a quantum state along phase
%
%   v=quadrature_variance(state, phase, S) calculates the x-quadrature variance of
%   the state along the phase.  state may be a vector or density matrix in the Fock basis.
%   This software assigns variance 1/2 to the vacuum state.  S is an
%   optional input containing the tables made by init_tables.  If S is not
%   present, it will be calculated.

if ~exist('S', 'var')
    photons=length(state)-1;
    S = init_tables(photons);
end

pe = phase_evolution(phase, S.photons);

if isvector(state)
    state = pe*state;
    expX = real(state'*S.expXMatrix*state);
    expXSquared = real(state'*S.expXSquaredMatrix*state);
else
    state = pe*state*pe';
    expX = real(trace_of_product(state, S.expXMatrix));
    expXSquared = real(trace_of_product(state, S.expXSquaredMatrix));
end

v = expXSquared-(expX)^2;

end

