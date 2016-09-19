function F = fidelity(rho, sigma)
% Calculates the fidelity of states rho and sigma, represented as density
% matrices. Similar to the formula from Nielsen and Chuang, p. 409, except
% squared.  This calculates the probability fidelity.
% 
% This can handle both density matrices and pure state column vectors.

if isvector(rho) && isvector(sigma)
    F = abs(rho'*sigma).^2;
elseif isvector(rho)
    F = rho'*sigma*rho;
elseif isvector(sigma);
    F = sigma'*rho*sigma;
else
    rho_half = sqrtm(rho);
    F = (trace(sqrtm(rho_half * sigma * rho_half)))^2;
end

if imag(F)>10^(-6)
    warning('Tomography:UnexpectedImaginary','fidelity has imaginary part larger than 10^-6')
end

F = real(F);

% TODO: when rho or sigma is rank 1, find its eigenvector and use that to
% compute fidelity.  Also perform sqrtm opration on matrix that has smaller
% condition number.