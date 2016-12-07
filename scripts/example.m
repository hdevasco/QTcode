clear all; 
clc;

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

% We choose a random list of phases at which our homodyne detecotr will
% measure this state.
nMeasurements = 20000;% Number of measurements
m = 20;% Number of equally spaced angles
                angles = pi*[0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19]/m; 
               angles = repmat(angles,1,ceil (nMeasurements/m)); 
                angles = angles(1:nMeasurements).';
                
% The homodyne detector has efficiency 90 %.
etaDetector = 0.9;
% Now, we make the measurements of state rho.  Notice we have to specify
% maximum and minimum possible measurement results -7 and 7.
samples = homodyne_samples(-7,7,etaDetector,angles,rho,S);

% In matrix H, each column equals the values of the quadrature measurements of the samples matrix for each angle separately.
H = zeros(1000,20);
% In matrix A, each angle is repeated a thousand times in each column in order to transform it into a column vector for the construction of M.
A = zeros(1000,20);
% In matrix N, we have the count number of the bins of the histograms of each equally spaced angle.
N = zeros(1000,20);
% In matrix C, each column equals the centers of the boxes of the histograms of each angle separately, both placed in 1000 bins.
C = zeros(1000,20);


    for (i=1:20)
        
        A(:,i) = samples((i:20:end),1);
        
        H(:,i) = samples((i:20:end),2); 
  % The histcounts function calculates the values at the edges and the number of counts in each bin for each angle at a time.
       
         [N(:,i),edges] = histcounts(H(:,i),1000);
  % "d" is the increment to be added to each edge to obtain the centers of the boxes
  
       d = diff(edges)/2;
  % Centers calculates the centers of the boxes in each on each histogram at a time
  
       centers = edges(1:end-1)+d;
       
        C(:,i)= centers';
        
        figure;

        title('h(i)');
     %"h (i)" plots the histograms of the quadrature measurements of each angle separately by placing them in 1000 bins.
     
        h(i) = histogram(H(:,i),1000)
        
    end
  % In matrix M, the columns are respectively (angles, centers of bins, number of counts in each bin)
  
   M=[A(:),C(:),N(:)];

%[homodyne_hist] = homodyne_hist(samples,bidWidth)
% Structure containing the POVM element corresponding to each measurement
% result.  Note that the POVMs are not pure projectors.  The homodyne
% detector's efficiency has been included in the computation of the POVMs.
%Povms = make_measurement_struct(samples,etaDetector,S);

% Now we will use the R*rho*R algorithm until we have done 2000 iterations
% or we reach stoppingCriterion 0.01 (whichever happens first).
% stoppingCriterion is an upper bound on the difference between the true
% maximum log-likelihood and the log-likelihood of that iterations's state.
%maxIterations = 2000;%20 * 1000 (The measurements of each angle are placed in 1000 bins)
%stoppingCriterion = 0.01;
%[rhoML1, Diagnostics] = rrhor_optimization(samples, S, etaDetector, 0, maxIterations, stoppingCriterion, []);
% output is the maximum likelihood state, rhoML, and a big structure
% Diagnostics containing information about each iterations's progress.

% Rather than using R*rho*R, we can use a combination of R*rho*R followed
% by iterations of our new algorithm, the regularized gradient ascent.
% This is usually faster, especially if you want to be very close to the
% true maximum likelihood state.
%[rhoML2, Diagnostics ] = combined_optimization( samples, S, etaDetector, 0, maxIterations, stoppingCriterion);

% Fidelity of true state and estimate
%fidelity(rhoML2, rho)


% We can plot the Wigner function of the the maximum likelihood state.
%wignerStepSize = 0.1;
%[x,p] = meshgrid(-6:wignerStepSize:6,-6:wignerStepSize:6);
%wignerML = wigner(rhoML2, x,p);
%figure;
%pcolor(x,p,wignerML)