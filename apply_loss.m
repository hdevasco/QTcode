function rhoEta = apply_loss( state, eta, S, opt )
% apply_loss transforms rho -> rhoEta after loss
%   input: state = state which will suffer loss.  It may be a pure column
%                  vector or a density matrix.  It may also be a POVM or a
%                  pure projector represented as a column vector.
%          eta = (scalar) transmissivity of medium through which rho passes
%                or eta may be a 3D array containing the elements of the 
%                loss operator sum representation. (Actually eta may be
%                the operator sum elements of any quantum operation).
%          S = init_tables(photons)
%          if opt = 'DoHermitianConjugate', it will perform hermitian
%          conjugate (actually just transpose because eta should be real) 
%          on the loss operator.  This is useful if the input eta
%          is the operator elements for application to a state, but you
%          want to apply the loss to an operator or POVM (or vice-verca).
%   output: rhoEta = new density matrix after loss

photons = S.photons;

% Decide whether we need to calculate E, and check for appropriate input
% for eta
if isscalar(eta)
    if isreal(eta) && 0<=eta && eta<=1
        E = loss_operation(eta, S);
    else
        error('Tomography:lossInvalid','efficiency must be a real number between 0 and 1')
    end
elseif ndims(eta) == 3 && size(eta,1) == size(eta,2) && size(eta,2) == size(eta,3)
    E = eta;
else
    error('Tomography:lossInvalid','loss must be a scalar or cube array')
end

% If the loss is applied to a POVM, we need to multiply the POVM by the
% hermitian conjugates of E(:,:,n).  E is real, so we only need to 
% transpose.
if nargin == 4 && strcmpi(opt, 'DoHermitianConjugate')
        E = permute(E, [2,1,3]);
end

% if the input state is a pure vector, the calculation can be much faster
if isvector(state)
    % The following four lines compute sum_n e_n|state><state|e_n^dagger.
    ETranspose = permute(E, [2, 1, 3]);
    EReshape = reshape(ETranspose, [S.dimHilbertSpace, (S.dimHilbertSpace).^2]);
    rhoEta = reshape(state'*EReshape,[S.dimHilbertSpace, S.dimHilbertSpace]);
    rhoEta = conj(rhoEta)*rhoEta.';
elseif ndims(state) == 2 && size(state,1) == size(state,2)
    % Make hermitian conjugate of each of e(:,:,n)
    % e is real, so I only need to transpose.
    eHC = permute(E,[2,1,3]);

    % The following three lines do matrix multiplication for each page of e, so
    % we get a 3D array each of whose pages is e(:,:,n)*rho*E_hc(:,:,n)
    % I had an other version of this function which accomplished the matrix
    % multiplication using a for loop, but it was slower than using arrayfun.
    % I also tried a user created toolbox called Multiprod, but it was MUCH
    % slower than both arrayfun and a for loop.
    rhoEta = arrayfun(@(n) E(:,:,n)*state*eHC(:,:,n), 1:photons+1,'uniformOutput',false);
    rhoEta = cell2mat(rhoEta);
    rhoEta = reshape(rhoEta,[photons+1,photons+1,photons+1]);

    %add the elements of the operator sum
    rhoEta = sum(rhoEta,3);
    
end

% TODO maybe this could be improved with bsxfun