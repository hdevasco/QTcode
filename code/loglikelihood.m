function ll = loglikelihood(rho, Measurements, tprl)
% loglikelihood computs the loglikelihood function for rho, measurements.
%
%   ll=loglikelihood(rho, Measurements) computes the logarithm of the
%   likelihood that measurements would be caused by density matrix rho.  Measurements is
%   the structre of quadrature measurement results and povms. 
%
%   ll = loglikelihood(rho, Measurements, tprl) uses the data in tprl to
%   speed the calculation.  tprl is the list of trace(rho*povm) made by
%   tr_povm_rho_list.

% if necessary make the trpl using the POVMlist
if nargin == 2
    tprl = tr_povm_rho_list(rho, Measurements);
end

ll = sum(Measurements.measurementArray(:,3).*log(tprl));

if imag(ll)/abs(ll)>10^(-6)
    warning('Tomography:UnexpectedImaginary','loglikelihood has imaginary part larger than 10^-6')
end

ll = real(ll);
