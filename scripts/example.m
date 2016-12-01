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
%Now let's form histograms with the quadrature values for each angle separately. 
%Each vector below represents the quadrature measurements for a given equally spaced angle of 0 to pi.

% Histogram for the quadrature measurements of the angle 0
a = samples(1:20:end,2);
figure;
title('h1');
[h1] = histogram(a,1000)
[N1,edges1] = histcounts(a,1000)
%[counts1,centers1] = hist(a)

% Histogram for the quadrature measurements of the angle 0.1571
b = samples(2:20:end,2);
figure;
title('h2');
[h2] = histogram(b,1000)
[N2,edges2] = histcounts(b,1000)
%[counts2,centers2] = hist(b)

% Histogram for the quadrature measurements of the angle 0.3142
c = samples(3:20:end,2);
figure;
title('h3');
[h3] = histogram(c,1000)
[N3,edges3] = histcounts(c,1000)
%[counts3,centers3] = hist(c)

% Histogram for the quadrature measurements of the angle 0.4712
d = samples(4:20:end,2);
figure;
title('h4');
[h4] = histogram(d,1000)
[N4,edges4] = histcounts(d,1000)
%[counts4,centers4] = hist(d)

% Histogram for the quadrature measurements of the angle 0.6283
e = samples(5:20:end,2);
figure;
title('h5');
[h5] = histogram(e,1000)
[N5,edges5] = histcounts(e,1000)
%[counts5,centers5] = hist(e)

% Histogram for the quadrature measurements of the angle 0.7854
f = samples(6:20:end,2);
figure;
title('h6');
[h6] = histogram(f,1000)
[N6,edges6] = histcounts(f,1000)
%[counts6,centers6] = hist(f)

% Histogram for the quadrature measurements of the angle 0.9425
g = samples(7:20:end,2);
figure;
title('h7');
[h7] = histogram(g,1000)
[N7,edges7] = histcounts(g,1000)
%[counts7,centers7] = hist(g)

% Histogram for the quadrature measurements of the angle 1.0996
h = samples(8:20:end,2);
figure;
title('h8');
[h8] = histogram(h,1000)
[N8,edges8] = histcounts(h,1000)
%[counts8,centers8] = hist(h)

% Histogram for the quadrature measurements of the angle 1.2566
i = samples(9:20:end,2);
figure;
title('h9');
[h9] = histogram(i,1000)
[N9,edges9] = histcounts(i,1000)
%[counts9,centers9] = hist(i)

% Histogram for the quadrature measurements of the angle 1.4137
j = samples(10:20:end,2);
figure;
title('h10');
[h10] = histogram(j,1000)
[N10,edges10] = histcounts(j,1000)
%[counts10,centers10] = hist(j)

% Histogram for the quadrature measurements of the angle 1.5708
k= samples(11:20:end,2);
figure;
title('h11');
[h11] = histogram(k,1000)
[N11,edges11] = histcounts(k,1000)
%[counts11,centers11] = hist(k)

% Histogram for the quadrature measurements of the angle 1.7279
l = samples(12:20:end,2);
figure;
title('h12');
[h12] = histogram(l,1000)
[N12,edges12] = histcounts(l,1000)
%[counts12,centers12] = hist(l)

% Histogram for the quadrature measurements of the angle 1.8850
m = samples(13:20:end,2);
figure;
title('h13');
[h13] = histogram(m,1000)
[N13,edges13] = histcounts(m,1000)
%[counts13,centers13] = hist(m)

% Histogram for the quadrature measurements of the angle 2.0420
n = samples(14:20:end,2);
figure;
title('h14');
[h14] = histogram(n,1000)
[N14,edges14] = histcounts(n,1000)
%[counts14,centers14] = hist(n)

% Histogram for the quadrature measurements of the angle 2.1991
o = samples(15:20:end,2);
figure;
title('h15');
[h15] = histogram(o,1000)
[N15,edges15] = histcounts(o,1000)
%[counts15,centers15] = hist(o)

% Histogram for the quadrature measurements of the angle 2.3562
p = samples(16:20:end,2);
figure;
title('h16');
[h16] = histogram(p,1000)
[N16,edges16] = histcounts(p,1000)
%[counts16,centers16] = hist(p)

% Histogram for the quadrature measurements of the angle 2.5133
q = samples(17:20:end,2);
figure;
title('h17');
[h17] = histogram(q,1000)
[N17,edges17] = histcounts(q,1000)
%[counts17,centers17] = hist(q)

% Histogram for the quadrature measurements of the angle 2.6704
r = samples(18:20:end,2);
figure;
title('h18');
[h18] = histogram(18,1000)
[N18,edges18] = histcounts(r,1000)
%[counts18,centers18] = hist(r)

% Histogram for the quadrature measurements of the angle 2.8274
t = samples(19:20:end,2);
figure;
title('h19');
[h19] = histogram(t,1000)
[N19,edges19] = histcounts(t,1000)
%[counts19,centers19] = hist(t)

% Histogram for the quadrature measurements of the angle 2.9845
u = samples(20:20:end,2);
figure;
title('h20');
[h20] = histogram(u,1000)
[N20,edges20] = histcounts(u,1000)
%[counts20,centers20] = hist(u)


%[homodyne_hist] = homodyne_hist(samples,bidWidth)
% Structure containing the POVM element corresponding to each measurement
% result.  Note that the POVMs are not pure projectors.  The homodyne
% detector's efficiency has been included in the computation of the POVMs.
Povms = make_measurement_struct(samples,etaDetector,S);

% Now we will use the R*rho*R algorithm until we have done 2000 iterations
% or we reach stoppingCriterion 0.01 (whichever happens first).
% stoppingCriterion is an upper bound on the difference between the true
% maximum log-likelihood and the log-likelihood of that iterations's state.
maxIterations = 2000;%20 * 1000 (The measurements of each angle are placed in 1000 bins)
stoppingCriterion = 0.01;
[rhoML1, Diagnostics] = rrhor_optimization(samples, S, etaDetector, 0, maxIterations, stoppingCriterion, []);
% output is the maximum likelihood state, rhoML, and a big structure
% Diagnostics containing information about each iterations's progress.

% Rather than using R*rho*R, we can use a combination of R*rho*R followed
% by iterations of our new algorithm, the regularized gradient ascent.
% This is usually faster, especially if you want to be very close to the
% true maximum likelihood state.
[rhoML2, Diagnostics ] = combined_optimization( samples, S, etaDetector, 0, maxIterations, stoppingCriterion);

% Fidelity of true state and estimate
fidelity(rhoML2, rho)


% We can plot the Wigner function of the the maximum likelihood state.
%wignerStepSize = 0.1;
%[x,p] = meshgrid(-6:wignerStepSize:6,-6:wignerStepSize:6);
%wignerML = wigner(rhoML2, x,p);
%figure;
%pcolor(x,p,wignerML)