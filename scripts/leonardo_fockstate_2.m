tic;
clc;
clear;
% This code simulates quantum state tomography using histograms in quadrature measurements.
% The histogram bin width is estimated for each simulation by the estimated average number of photons. 
% In the end, the fidelity between the true state and the rebuilt state is calculated.
% Fock State for n=[4,6,8,10]

%Maximum number of photons

maxPhotonNumber = [15];

% Number of measurements
numMeasurements    = 20000;


    
% Number of simulations
numSim             = 100;

etaDetector        = 0.9;
maxIterations      = 2000;
stoppingCriterion  = 0.01;
 n                = [4,6,8,10];
etaState           = 0.95; 

% Number of angles equally spaced from 0 to pi
numAngles          = 20;

                
 for i=1:length(maxPhotonNumber),
    for k=1:length(numMeasurements),
         for j=1:length(n),
            
            fileName1 = ['maxPhNum',num2str(maxPhotonNumber(i)),'n',num2str(n(j)),'nM',num2str(numMeasurements(k)),'.mat'];
            
            if exist(fileName1,'file') == 2,
                load(fileName1);

            else
                fML2                = zeros(numSim,1);
                fScott              = zeros(numSim,1);
                fHistogram          = zeros(numSim,1);
                fHistogramPsi       = zeros(numSim,1);
                fML2Psi             = zeros(numSim,1);
                fScottPsi           = zeros(numSim,1);
                fidelityDiff        = zeros(numSim,1);
                fidelityDiff2       = zeros(numSim,1);
                 n_bar              = zeros(numSim,1);
                 deltaq             = zeros(numSim,1);
                 n_RhoHistogram     = zeros(numSim,1);
                 n_RhoScott         = zeros(numSim,1);
                 n_RhoML2           = zeros(numSim,1);
                 t                  = 1;
                
                save(fileName1);
            end
            
            while t <= numSim,
                fprintf(['Iteraction ',num2str(t), '\n']);
                
              angles = pi*(0:numAngles-1)/numAngles;
              angles = repmat(angles,1, ceil(numMeasurements/numAngles));
              angles = angles(1:numMeasurements)';
              %generate state 

               tic;
                S   = init_tables(maxPhotonNumber(i));
                psi = generate_fock_vector(n(j), maxPhotonNumber(i));
                Rho = apply_loss(psi,etaState,S);
            
               Samples = homodyne_samples(-7,7,etaDetector,angles,Rho,S);
                
                % average number of photons estimated ( \overline{\langle
                % \hat{n} \rangle}) \overline{\langle \hat{n} \rangle}
                n_bar(t) = n_quadrature(Samples,numMeasurements);
               
                % We determined the histogram box width from the estimated average number of photons
                
                 deltaq(t) = pi/(2*sqrt(2*n_bar(t)+1));

               
                [RhoML2, Diagnostics] = combined_optimization( Samples, S, etaDetector, 0, maxIterations, stoppingCriterion);
                
                % estimates the average number of photons of the
                % reconstructed state RhoML2
                 n_RhoML2(t) = mean_photons(RhoML2);
                
              % Calculates the fidelity of RhoML2 and Rho
                fML2(t) = fidelity(RhoML2, Rho);
                    
                % Constructs the MHistogram matrix using the given number of bins (case option> 6 in the matrix_histogram function)
               
 
                %  Constructs RhoHistogram using the combined_optimization
                 MHistogram = matrix_histogram(Samples, 7,'integral', deltaq(t));
                
                 [RhoHistogram, Diagnostics] = combined_optimization( MHistogram, S, etaDetector, 0, maxIterations, stoppingCriterion);
                 
              
                % estimates the average number of photons of the
               % reconstructed state RhoHistogram
                 n_RhoHistogram(t) = mean_photons(RhoHistogram);
                 
                % Constructs the MScott array using the optimal width method (case option = 8 in the matrix_histogram function)
                tic;
                MScott = matrix_histogram(Samples, 8,'integral');
   
                % Constructs RhoScott using the combined_optimization
                [RhoScott, Diagnostics] = combined_optimization( MScott, S, etaDetector, 0, maxIterations, stoppingCriterion);
               
                % estimates the average number of photons of the
                % reconstructed state RhoScott
                 n_RhoScott(t) = mean_photons(RhoScott);
                
                % Calculates the fidelity of RhoHistogram and Rho
                fHistogram(t)       = fidelity(RhoHistogram, Rho);
            
                % Calculates the fidelity of RhoSCott and Rho
                fScott(t)           = fidelity(RhoScott, Rho);     
                
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
                meanFScott          = mean(fScott);
                meanFML2Psi         = mean(fML2Psi);
                meanFHistogramPsi   = mean(fHistogramPsi);
                meanFScottPsi         = mean(fScottPsi); 
                mean_n_bar              = mean(n_bar); % average number of estimated photons
                mean_deltaq             = mean(deltaq);
                mean_photon_RhoML2      = mean(n_RhoML2);
                mean_photon_RhoHistogram = mean(n_RhoHistogram);
                mean_photon_RhoScott    = mean(n_RhoScott);
                std_n_RhoML2            = std(n_RhoML2);
                std_n_RhoHistogram      = std(n_RhoHistogram);
                std_n_RhoScott          = std(n_RhoScott);
                stdFHistogram           = std(fHistogram);        % Standard deviation of fHistogram
                stdFML2                 = std(fML2);              % Standard deviation of fML2
                stdFScott               = std(fScott);            % Standard deviation of fScott
                std_n_bar               = std(n_bar);
                std_deltaq              = std(deltaq);
               
                home;
                fprintf('>> Progress: %.2f%%\n', t/numSim*100);
                t=t+1;
                
                save(fileName1);
            end
        end
    end
end   