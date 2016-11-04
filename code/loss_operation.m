function E = loss_operation(eta, S, opt)
% returns loss operation elements in the operator-sum represtnation.  These
% are used to calculate new density matrix after light travels through
% medium with transmissivity eta.  rho_new = Sum_k^{maxphoton} E(n) rho
% E(n)'.  E is a three dimensional array, and E(:,:,n) is the nth operator
% in the sum. S is the data structre with things like binomial coefficients
% and Hermite polynomials.  If you intent to apply E to an operator or
% POVM, then give opt = 'ToPOVM', the loss_operation will return the
% hermitian conjugate of the operator elements.
%
% to write out an expression for E(n) in latex:
% E(j,k,n) = sum_{l=0}^k \sqrt{binom{k}{l} \sqrt{\eta}^l \sqrt{1-\eta}^k-l
% <j|l><n|k-l>

E = zeros(S.dimHilbertSpace,S.dimHilbertSpace,S.dimHilbertSpace);

if nargin == 2
    % E(:,:,n) is a diagonal matrix, so I use diag to form it.  v is the
    % list of values that occupy the jth diagonal.
    for j = 0:S.photons
        v = sqrt(diag(S.binom,-j).'.*eta.^(0:S.photons-j).*(1-eta).^j);
        E(:,:,j+1) = diag(v,j);
    end
elseif nargin == 3 && strcmpi(opt, 'ToPOVM')
    for j = 0:S.photons
        v = sqrt(diag(S.binom,-j).'.*eta.^(0:S.photons-j).*(1-eta).^j);
        % Using -j in the diag function creates the transpose of diag(v,j).
        E(:,:,j+1) = diag(v,-j);
    end
end

%{
Another version I tried used arrayfun, but it was slightly slower
E = arrayfun(@(n) diag(v(n,S.photon,S.binom,eta),n),0:S.maxPhoton,'UniformOutput',false);
E = cell2mat(E);
E = reshape(E,[S.photon+1,S.photon+1,S.photon+1]);
%}

%{
An earlier version that did not use diag was much slower.  For reference, 
here it is

% 3D arrays that contain indexs to column, row, and depth
% k labels the rows, j labels the columns, and n is the depth
index_k = repmat(0:max_photon,[max_photon+1,1,max_photon+1]);
index_j = permute(index_k,[2,1,3]);
index_n = permute(index_k,[1,3,2]);

E = zeros(max_photon+1,max_photon+1,max_photon+1);

% Most elements of E will be zero, so we will only compute the nonzero
% elements, using logical indexing.  H contains 1's at each location where
% E may be nonzero.
H = index_j<=index_k & index_k-index_j==index_n;

eta_powers = zeros(max_photon+1,max_photon+1,max_photon+1);
one_minus_eta_powers = zeros(max_photon+1,max_photon+1,max_photon+1);
binomials = repmat(S.binom.',[1,1,max_photon+1]);

% We will replace the zeros with nonzero factors that we will later
% multiply to get E.
eta_powers(H) = eta.^index_j(H);
one_minus_eta_powers(H) = (1-eta).^(index_k(H)-index_j(H));

E(H) = sqrt(binomials(H).*eta_powers(H).*one_minus_eta_powers(H));
%}

