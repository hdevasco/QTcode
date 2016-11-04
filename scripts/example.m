% Here is an example using the tomography code to simulate measurements
% made on an optical "Schordinger cat state" and to find the maximum
% likelihood state for that set of measurements.

% The infinite dimensional state space for the harmonic oscillator will be
% represented in the photon number basis. 
% We will truncate the Hilbert space at maxPhotonNumber photons.
maxPhotonNumber = 10;

% First, pre-compute a lot of numbers, such as coefficients for Hermite
% polynomials, factorials, binomial coefficients.
S = init_tables(maxPhotonNumber);

% Make state vector for Schrodinger cat state.
alpha = 1;  % amplitude of coherent states in the superposition
phase = 0;  % phase between superposition
psi = generate_cat_vector(alpha, phase, S);
% The Schrodinger cat state suffers from some loss by passing through a
% medium with 80 % efficiency.
etaState = 0.8;
rho = apply_loss(psi,etaState,S);
% Now it must be represented by a density matrix, rho.

% We choose a random list of phases at which our homodyne detector will
% measure this state.
nMeasurements = 40000;
angles = 2*pi*rand(nMeasurements,1);
% The homodyne detector has efficiency 90 %.
etaDetector = 0.9;
% Now, we make the measurements of state rho.  Notice we have to specify
% maximum and minimum possible measurement results -7 and 7.
samples = homodyne_samples(-7,7,etaDetector,angles,rho,S);
% samples is a 3 x length(angles) array, each row containing (angle,
% quadrature, 1).  The "1" means that that angle, quadrature combination
% was observed once.  samples does not bin or histogram the data.  If you
% want to bin your angles and / or quadratures, you can prepare your own 3
% x N array with the number of times each angle / quadrature was observed
% in the third column of the array.

% Structure containing the POVM element corresponding to each measurement
% result.  Note that the POVMs are not pure projectors.  The homodyne
% detector's efficiency has been included in the computation of the POVMs.
Povms = make_measurement_struct(samples,etaDetector,S);

% Now we will use the R*rho*R algorithm until we have done 2000 iterations
% or we reach stoppingCriterion 0.01 (whichever happens first).
% stoppingCriterion is an upper bound on the difference between the true
% maximum log-likelihood and the log-likelihood of that iterations's state.
maxIterations = 2000;
stoppingCriterion = 0.2;
[rhoML1, Diagnostics1] = rrhor_optimization(samples, S, etaDetector, 0, maxIterations, stoppingCriterion, []);
% output is the maximum likelihood state, rhoML1, and a big structure
% Diagnostics1 containing information about each iterations's progress.

% Rather than using R*rho*R, we can use a combination of R*rho*R followed
% by iterations of our new algorithm, the regularized gradient ascent.
% This is usually faster, especially if you want to be very close to the
% true maximum likelihood state.
[rhoML2, Diagnostics2 ] = combined_optimization( samples, S, etaDetector, 0, maxIterations, stoppingCriterion);

% Fidelity of true state and estimate
fidelity(rhoML2, rho)


% We can plot the Wigner function of the the maximum likelihood state.
wignerStepSize = 0.1;
[x,p] = meshgrid(-6:wignerStepSize:6,-6:wignerStepSize:6);
wignerML = wigner(rhoML2, x,p);
pcolor(x,p,wignerML)