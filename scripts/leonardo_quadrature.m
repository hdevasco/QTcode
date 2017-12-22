tic;
clc;
clear;
% This code simulates quantum state tomography using histograms in quadrature measurements.
% Calculates the estimate for the mean number of photons using the quadrature measurements and the correction 
% by dividing this estimated value by the efficiency of the homodyne detector.
% We also calculate the mean of the mean numbers of the reconstructed density matrices using the histogram.

%Maximum number of photons

maxPhotonNumber = [10];

% Number of measurements
numMeasurements    = 20000;

    
% Number of simulations
numSim             = 100;

etaDetector        = 0.9;
maxIterations      = 2000;
stoppingCriterion  = 0.01;
alpha              = [1] ;
phase              = 0;
etaState           = 0.8;

% Number of angles equally spaced from 0 to pi
numAngles          = 20;

                
 for i=1:length(maxPhotonNumber),
    for k=1:length(alpha),
         for j=1:length(numMeasurements),
            
            fileName1 = ['maxPhNum',num2str(maxPhotonNumber(i)),'nM',num2str(numMeasurements(j)),'alpha',num2str(alpha(k)),'.mat'];
            
            if exist(fileName1,'file') == 2,
                load(fileName1);
            else
                fML2                = zeros(numSim,1);
                fScott              = zeros(numSim,1);
                fHistogram          = zeros(numSim,1);
                fHistogramC          = zeros(numSim,1);
                fHistogramPsi       = zeros(numSim,1);
                fML2Psi             = zeros(numSim,1);
                fScottPsi           = zeros(numSim,1);
                fidelityDiff        = zeros(numSim,1);
                fidelityDiff2       = zeros(numSim,1);
                timeML2             = zeros(numSim,1);
                n_bar               = zeros(numSim,1);
                n                   = zeros(numSim,1);
                deltaq              = zeros(numSim,1);
                deltaqC              = zeros(numSim,1);
                n_RhoHistogram      = zeros(numSim,1);
                n_RhoHistogramC      = zeros(numSim,1);
                n_RhoScott          = zeros(numSim,1);
                n_RhoML2            = zeros(numSim,1);
                 deltaqC              = zeros(numSim,1);
                n_RhoHistogramC      = zeros(numSim,1);
                n_RhoML2            = zeros(numSim,1);
                t                   = 1;
                
                save(fileName1);
            end
            
            while t <= numSim,
                fprintf(['Iteraction ',num2str(t), '\n']);
             
             % Divide the semicircle of phase space into equally spaced d-phases
              angles = pi*(0:numAngles-1)/numAngles;
              angles = repmat(angles,1, ceil(numMeasurements/numAngles));
              angles = angles(1:numMeasurements)';
              
               % generate state 
                S   = init_tables(maxPhotonNumber(i));
                 psi = generate_cat_vector(alpha(k), phase, S);
                 Rho = apply_loss(psi,etaState,S);
             
                 % Generates phase quadrature measurements
                Samples = homodyne_samples(-7,7,etaDetector,angles,Rho,S);

                
                % estimates the average number of photons using the quadrature measurements
                n(t) = n_quadrature(Samples,numMeasurements);
                
                % corrects the estimation for the average number of photons by dividing by the efficiency of the detector
                n_bar(t) = n(t)/etaDetector;
               
                % we determined the width of the histogram box using the estimation of the average number of photons
                deltaq(t) = pi/(2*sqrt(2*n(t)+1));
                
                % we determined the width of the histogram box using the estimation of the mean number of photons after correction
                deltaqC(t) = pi/(2*sqrt(2*n_bar(t)+1));
                
                
                [RhoML2, Diagnostics] = combined_optimization( Samples, S, etaDetector, 0, maxIterations, stoppingCriterion);
                
                
                % estimates the average number of photons of the
                % reconstructed state RhoML2
                n_RhoML2(t) = mean_photons(RhoML2);
                
              % Calculates the fidelity of RhoML2 and Rho
                fML2(t) = fidelity(RhoML2, Rho);
                    
                % Constructs the MHistogram matrix using the given number of bins (case option> 6 in the matrix_histogram function)
               
                MHistogram = matrix_histogram(Samples, 7,'integral', deltaq(t));
                
                % Constructs the MHistogramC matrix(with the width after correction in the average number of photons)
                MHistogramC = matrix_histogram(Samples, 7,'integral', deltaqC(t));
                
                %  Constructs RhoHistogram using the combined_optimization
                [RhoHistogram, Diagnostics] = combined_optimization( MHistogram, S, etaDetector, 0, maxIterations, stoppingCriterion);
               
                
                 %  Constructs RhoHistogram using the combined_optimization(with the width after correction in the average number of photons)
                [RhoHistogramC, Diagnostics] = combined_optimization( MHistogramC, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeRhoHistogram(t) = toc;
                
                % estimates the average number of photons of the
                % reconstructed state RhoHistogram
                n_RhoHistogram(t) = mean_photons(RhoHistogram);
                
                %(after correction in the average number of photons)
                n_RhoHistogramC(t) = mean_photons(RhoHistogramC);
                
                % Constructs the MScott array using the optimal width method (case option = 2 in the matrix_histogram function)
                
                MScott = matrix_histogram(Samples, 2,'integral');
   
                % Constructs RhoScott using the combined_optimization
                [RhoScott, Diagnostics] = combined_optimization( MScott, S, etaDetector, 0, maxIterations, stoppingCriterion);
               
                % estimates the average number of photons of the
                % reconstructed state RhoScott
                n_RhoScott(t) = mean_photons(RhoScott);
                
                fHistogram(t)       = fidelity(RhoHistogram, Rho);
                fHistogramC(t)       = fidelity(RhoHistogramC, Rho); 
                % Calculates the fidelity of RhoHistogram and Rho
                fScott(t)           = fidelity(RhoScott, Rho);     % Calculates the fidelity of RhoScott and Rho
               % fHistogramPsi(t)    = fidelity(RhoHistogram, psi); % Calculates the fidelity of RhoHistogram and Psi
               % fML2Psi(t)          = fidelity(RhoML2, psi);       % Calculates the fidelity of RhoML2 and Psi
                %fScottPsi(t)        = fidelity(RhoScott, psi);     % Calculates the fidelity of RhoScott and Rho
                
                fidelityDiff(t)     = fML2(t)-fHistogram(t);       % Difference between fidelities [fidelity(RhoML2,Rho)-fidelity(RhoHistogram,Rho)]
                fidelityDiff2(t)    = fML2(t)-fScott(t);           % Difference between fidelities [fidelity(RhoML2,Rho)-fidelity(RhoScott,Rho)]
                wHistogram          = mean(fidelityDiff);
                wHistogram2         = mean(fidelityDiff2);
                Dhistogram          = std(fidelityDiff);           % Standard deviation of the difference between the fidelities(fML2,fHistogram)
                Dhistogram2          = std(fidelityDiff2);         % Standard deviation of the difference between the fidelities(fML2,fScott)
                meanFML2            = mean(fML2);                  
                meanFHistogram      = mean(fHistogram);
                meanFHistogramC      = mean(fHistogramC);
                meanFScott          = mean(fScott);
                meanFML2Psi         = mean(fML2Psi);
                meanFHistogramPsi   = mean(fHistogramPsi);
                meanFScottPsi           = mean(fScottPsi);     
                meanTimeRhoML2          = mean(timeML2);
                mean_n_bar              = mean(n_bar);
                mean_n                  = mean(n);
                mean_deltaq             = mean(deltaq);
                mean_deltaqC             = mean(deltaqC);
                mean_photon_RhoML2      = mean(n_RhoML2);
                mean_photon_RhoHistogram = mean(n_RhoHistogram);
                mean_photon_RhoHistogramC = mean(n_RhoHistogramC);
                mean_photon_RhoScott    = mean(n_RhoScott);
                std_n_RhoML2            = std(n_RhoML2);
                std_n_RhoHistogram      = std(n_RhoHistogram);
                std_n_RhoHistogramC      = std(n_RhoHistogramC);
                std_n_RhoScott          = std(n_RhoScott);
                stdFHistogram           = std(fHistogram);        % Standard deviation of fHistogram
                stdFML2                 = std(fML2);              % Standard deviation of fML2
                stdFScott               = std(fScott);            % Standard deviation of fScott
                std_n_bar               = std(n_bar);
                std_n                   = std(n);
                std_deltaq              = std(deltaq);
                std_deltaqC              = std(deltaqC);
                home;
                fprintf('>> Progress: %.2f%%\n', t/numSim*100);
                t=t+1;
                
                save(fileName1);
            end
        end
    end
end   