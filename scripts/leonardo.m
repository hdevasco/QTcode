tic;
clc;
clear;
%Maximum number of photons
maxPhotonNumber    = 10;

% Number of measurements
numMeasurements    = 20000;

% Simulates the reconstruction of the quantum state using the Scott's Method and the bins number method determined for the histogram
numBins            = 100;
method             = 2;
option             = [method, numBins];

% Number of simulations
numSim             = 100;

etaDetector        = 0.9;
maxIterations      = 2000;
stoppingCriterion  = 0.01;
alpha              = 1;
phase              = 0;
etaState           = 0.8;

numAngles          = 20;

for i=1:length(maxPhotonNumber),
    for k=1:length(numMeasurements),
         for j=1:length(option),
            
            fileName1 = ['maxPhNum',num2str(maxPhotonNumber(i)),'op',num2str(option(j)),'nM',num2str(numMeasurements(k)),'.mat'];
            
            if exist(fileName1,'file') == 2,
                load(fileName1);
            else
                fML2                = zeros(numSim,1);
                fHistogram          = zeros(numSim,1);
                fHistogramPsi       = zeros(numSim,1);
                fML2Psi             = zeros(numSim,1);
                fidelityDiff        = zeros(numSim,1);
                timeML2             = zeros(numSim,1);
                timeMhistogram      = zeros(numSim,1);
                timeRhoHistogram    = zeros(numSim,1);
                t                   = 1;
                
                save(fileName1);
            end
            
            angles = pi*(0:numAngles-1)/numAngles;
            angles = repmat(angles,1, ceil(numMeasurements(k)/numAngles));
            angles = angles(1:numMeasurements(k))';
            
            while t <= numSim,
                fprintf(['Iteraction ',num2str(t), '\n']);
                
                tic;
                S   = init_tables(maxPhotonNumber(i));
                psi = generate_cat_vector(alpha, phase, S);
                Rho = apply_loss(psi,etaState,S);
                
                Samples = homodyne_samples(-7,7,etaDetector,angles,Rho,S);
                
                tic;
                [RhoML2, Diagnostics] = combined_optimization( Samples, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeML2(t) = toc;
                
                fML2(t) = fidelity(RhoML2, Rho);

                tic;
                MHistogram = matrix_histogram(Samples, option(j));
                timeMhistogram(t) = toc;
                
                tic;
                [RhoHistogram, Diagnostics] = combined_optimization( MHistogram, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeRhoHistogram(t) = toc;
                
                fHistogram(t)       = fidelity(RhoHistogram, Rho);
                fHistogramPsi(t)    = fidelity(RhoHistogram, psi);
                fML2Psi(t)          = fidelity(RhoML2,psi);
                
                fidelityDiff(t)     = fML2(t)-fHistogram(t);
                wHistogram          = mean(fidelityDiff);
                Dhistogram          = std(fidelityDiff);
                meanfML2                = mean(fML2);
                meanFHistogram      = mean(fHistogram);
                meanFML2psi         = mean(fML2Psi);
                meanFHistogramPsi       = mean(fHistogramPsi);
                meanTimeML2             = mean(timeML2);
                meanTimeMHistogram      = mean(timeMhistogram);
                meanTimeRhoHistogram    = mean(timeRhoHistogram);
                stdFHistogram           = std(fHistogram);
                stdFML2                 = std(fML2);
                
                home;
                fprintf('>> Progress: %.2f%%\n', t/numSim*100);
                t=t+1;
                
                save(fileName1);
            end
        end
    end
end