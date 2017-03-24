tic;
clc;
clear;
%Maximum number of photons
maxPhotonNumber    = 10;

% Number of measurements
numMeasurements    = 20000;

% Simulates the reconstruction of the quantum state using the Scott's Method and the bins number method determined for the histogram
numBins            = 100;
option             = [numBins];

% Number of simulations
numSim             = 100;

etaDetector        = 0.9;
maxIterations      = 2000;
stoppingCriterion  = 0.01;
alpha              = 1;
phase              = 0;
etaState           = 1/2; % antes era 0.8

% Number of angles equally spaced from 0 to pi
numAngles          = 20;

for i=1:length(maxPhotonNumber),
    for k=1:length(numMeasurements),
         for j=1:length(option),
            
            fileName1 = ['maxPhNum',num2str(maxPhotonNumber(i)),'op',num2str(option(j)),'nM',num2str(numMeasurements(k)),'.mat'];
            
            if exist(fileName1,'file') == 2,
                load(fileName1);
            else
                fML2                = zeros(numSim,1);
                fScott              = zeros(numSim,1);
                fHistogram          = zeros(numSim,1);
                fHistogramVacumm    = zeros(numSim,1);
                fML2Vacumm          = zeros(numSim,1);
                fScottVacumm        = zeros(numSim,1);
                fidelityDiff        = zeros(numSim,1);
                fidelityDiff2       = zeros(numSim,1);
                timeML2             = zeros(numSim,1);
                timeMhistogram      = zeros(numSim,1);
                timeMScott          = zeros(numSim,1);
                timeRhoHistogram    = zeros(numSim,1);
                timeRhoScott        = zeros(numSim,1);
                t                   = 1;
                ratioSwitch     = 'true variance';
                save(fileName1);
            end
            
            angles = pi*(0:numAngles-1)/numAngles;
            angles = repmat(angles,1, ceil(numMeasurements(k)/numAngles));
            angles = angles(1:numMeasurements(k))';
            
            while t <= numSim,
                fprintf(['Iteraction ',num2str(t), '\n']);
                
                tic;
                S   = init_tables(maxPhotonNumber(i));
                
                % novo estado pedido Dr.Glancy, estado de vácuo comprimido
                % com variância de 3/4 do vácuo
                varianceOrRatio = 3/4;
                state = generate_squeezed_vacuum_vector(varianceOrRatio, maxPhotonNumber, ratioSwitch);
                
                Rho = apply_loss(state,etaState,S);
                
                Samples = homodyne_samples(-7,7,etaDetector,angles,Rho,S);
                
                % Constructs RhoML2 using the combined_optimization
                tic;
                [RhoML2, Diagnostics] = combined_optimization( Samples, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeML2(t) = toc;
                
                % Calculates the fidelity of RhoML2 and Rho
                fML2(t) = fidelity(RhoML2, Rho);
                    
                % Constructs the MHistogram matrix using the given number of bins (case option> 6 in the matrix_histogram function)
                tic;
                MHistogram = matrix_histogram(Samples, option(j));
                timeMhistogram(t) = toc;
                
                % Constructs the MScott array using the optimal width method (case option = 2 in the matrix_histogram function)
                tic;
                MScott = matrix_histogram(Samples, 2);
                timeMScott(t) = toc;
                
                %  Constructs RhoHistogram using the combined_optimization
                tic;
                [RhoHistogram, Diagnostics] = combined_optimization( MHistogram, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeRhoHistogram(t) = toc;
                
                % Constructs RhoScott using the combined_optimization
                tic;
                [RhoScott, Diagnostics] = combined_optimization( MScott, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeRhoScott(t) = toc;
                
                fHistogram(t)       = fidelity(RhoHistogram, Rho); % Calculates the fidelity of RhoHistogram and Rho
                fScott(t)           = fidelity(RhoScott, Rho);     % Calculates the fidelity of RhoScott and Rho
                fHistogramVacumm(t)    = fidelity(RhoHistogram, state); % Calculates the fidelity of RhoHistogram and Psi
                fML2Vacumm(t)          = fidelity(RhoML2, state);       % Calculates the fidelity of RhoML2 and Psi
                fScottVacumm(t)        = fidelity(RhoScott, state);     % Calculates the fidelity of RhoScott and Rho
                
                fidelityDiff(t)     = fML2(t)-fHistogram(t);       % Difference between fidelities [fidelity(RhoML2,Rho)-fidelity(RhoHistogram,Rho)]
                fidelityDiff2(t)    = fML2(t)-fScott(t);           % Difference between fidelities [fidelity(RhoML2,Rho)-fidelity(RhoScott,Rho)]
                wHistogram          = mean(fidelityDiff);
                wHistogram2         = mean(fidelityDiff2);
                Dhistogram          = std(fidelityDiff);           % Standard deviation of the difference between the fidelities(fML2,fHistogram)
                Dhistogram2         = std(fidelityDiff2);         % Standard deviation of the difference between the fidelities(fML2,fScott)
                meanFML2            = mean(fML2);                  
                meanFHistogram      = mean(fHistogram);
                meanFScott          = mean(fScott);
                meanFML2Vacumm         = mean(fML2Vacumm);
                meanFHistogramVacumm    = mean(fHistogramVacumm);
                meanFScottVacumm        = mean(fScottVacumm);     
                meanTimeRhoML2          = mean(timeML2);
                meanTimeMHistogram      = mean(timeMhistogram);
                meanTimeMScott          = mean(timeMScott);
                meanTimeRhoHistogram    = mean(timeRhoHistogram);
                meanTimeRhoScott        = mean(timeRhoScott);
                stdFHistogram           = std(fHistogram);        % Standard deviation of fHistogram
                stdFML2                 = std(fML2);              % Standard deviation of fML2
                stdFScott               = std(fScott);            % Standard deviation of fScott
                home;
                fprintf('>> Progress: %.2f%%\n', t/numSim*100);
                t=t+1;
                
                save(fileName1);
            end
        end
    end
end