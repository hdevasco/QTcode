function catPsi = generate_cat_vector(alpha, theta, maxPhotons)
% Creates Schrodinger cat state
%    catPsi = generate_cat_vector(alpha, theta, maxPhoton) returns the pure
%    state vector for a Schrodinger cat state (namely, 
%    e^{i*theta}|-alpha>+|alpha>) in the photon number basis. maxPhoton is
%    the photon number at which the Hilbert space is truncated or the table
%    made by init_tables.
% 

if isstruct(maxPhotons)
    maxPhotons = maxPhotons.photons;
end

catPsi = (exp(1i .* theta) .* (-alpha).^(0:maxPhotons) + alpha.^(0:maxPhotons)) ./ sqrt(factorial(0:maxPhotons));
catPsi = catPsi.';

normalization = exp(-abs(alpha).^2./2)./sqrt(2.*(1+exp(-2*abs(alpha)^2).*cos(theta)));
catPsi = normalization .* catPsi;

catPsi = normalize(catPsi,'check');