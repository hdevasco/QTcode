function R = make_r_struct(rho, Measurements, linearForm)
%makes structure containing data related to current density matrix
%   make_r_struct(rho, Measurements) makes a structrue containing data
%   related to the current density matrix rho.  It requires the structure
%   of Measurement data, created by make_measurement_struct.
%
%   make_r_struct(rho, Measurements, linearForm) includes the linearForm
%   for optimizing L(rho) + trace(rho*linearForm).  linearForm must be a
%   matrix with the same size as rho.


    R.rho = rho;
    R.tprl = tr_povm_rho_list(rho, Measurements);
    R.r = r_or_gradient(R.rho, Measurements, R.tprl, 'R');
    R.loglike = loglikelihood(R.rho, Measurements, R.tprl);

    if exist('linearForm', 'var')
        R.r = R.r + linearForm;
        if ismember(inf,R.r) 
            R.go = false; 
            R.stop=0;
        else
            R.stop = real(max_eig(R.r) - Measurements.nTotalMeasurements - trace_of_product(linearForm, rho));
            % note Measurements.nTotalMeasurements == trace(r*rho)
        end
            R.linearTerm = real(trace_of_product(rho, linearForm));
            R.objective = R.loglike + R.linearTerm;

    else
        R.stop = real(max_eig(R.r)) - Measurements.nTotalMeasurements;
    end

    end
