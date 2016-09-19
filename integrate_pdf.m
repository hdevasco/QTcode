function prob = integrate_pdf(rho, photons, eta, phase, limit, n)

S = init_tables(photons);
stepSize = 2*limit/n;
measurementArray = [repmat(phase,n+1,1),(-limit:stepSize:limit)'];
Measurements = make_measurement_struct(measurementArray, eta, S);
tprl = tr_povm_rho_list(rho, Measurements);
steps = size(tprl,1);

prob = sum(tprl(2:(steps-1)));
prob = prob + 0.5*(tprl(1)+tprl(end));
prob = stepSize*prob;

end

