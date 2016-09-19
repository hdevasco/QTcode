function tprl = tr_povm_rho_list( rho, Measurements)
%makes a list of trace(povm*rho) for every POVM in povmArray
%   tr_povm_rho_list(rho, Measurements) makes a list of trace(povm*rho)
%   for every povm in povmArray
%   input: rho = n-by-n density matrix;
%          Measurements = structre containing measurement results
%          and povms, made by make_measurement_struct
%   output: tprl = number_of_measurements by 1 list of trace(rho*povm) for
%                  each different measurement


% for n = 1:Measurements.nMeasurements
%     tprl(n) = trace_of_product(Measurements.povmArray(:,:,n),rho);
% end

tprl = bsxfun(@times, Measurements.povmArray, rho.');
tprl = sum(sum(tprl));
tprl = tprl(:);

end
