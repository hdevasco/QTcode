function rG = r_or_gradient( rho, Measurements, tprl, option)
%R matrix for R*rho*R iterations or the gradient of the log-likelihood function
%   r_or_gradient(rho, Measurements, tprl, option) is the
%   derivatives of the log-likelihood function with respect to the
%   elements of the density matrix.
%   inputs: rho = density matrix
%           Measurements = structre of measurements and POVMs
%           tprl = optional list of trace(povm*rho) created by 
%                  tr_povm_rho_list
%           option = (mandatory) 'R' if we are calculating the R
%                    matrix, 'gradient' if we want the gradient of
%                    the log-likelihood function
%   output: rG = gradient matrix. g(a,b) = derivative of the log-likelihood
%               function with respect to rho(a,b).

if nargin == 3
    option = tprl;
    tprl = tr_povm_rho_list(rho, Measurements);
end

itprl = Measurements.measurementArray(:,3)./tprl; % array of n/Tr(rho*P), where n is number of measurements for measurement operator P

% rG = zeros(size(Measurements.povmArray));
% for n = 1:Measurements.nMeasurements
%     rG(:,:,n) = Measurements.povmArray(:,:,n).*tprl(n);
% end

itprl = reshape(itprl, [1, 1, Measurements.nMeasurements]);
rG = bsxfun(@times, itprl, Measurements.povmArray);

rG = sum(rG, 3);

if strcmp(option,'gradient')
    rG = rG.';
end

end