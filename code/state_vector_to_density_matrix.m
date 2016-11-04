function rho  = state_vector_to_density_matrix( psi )
%STATE_VECTOR_TO_DENSITY_MATRIX transforms pure state vector to density
%matrix
% input: psi = pure state column vector in Fock basis
% output: rho = density matrix in Fock basis

rho = psi*psi';
