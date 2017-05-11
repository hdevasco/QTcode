tic;
clc;
clear;

%Maximum number of photons
Mph = [10];


% Specified Bin Width
deltaq =[0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.2, 0.34]
numSim             = 100;
numMeasurements = 20000;
etaDetector        = 0.9;
maxIterations      = 2000;
stoppingCriterion  = 0.01;
alpha              = 1;
phase              = 0;
etaState           = 0.8;

% Number of angles equally spaced from 0 to pi
numAngles          = 20;

for i=1:length(Mph),
    for k=1:length(numMeasurements),
         for j=1:length(deltaq),
            
            fileName1 = ['maxPhNum',num2str(Mph(i)),'BinWidth',num2str(deltaq(j)),'nM',num2str(numMeasurements(k)),'.mat'];
            
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
                timeML2             = zeros(numSim,1);
                timeMhistogram      = zeros(numSim,1);
                timeMScott          = zeros(numSim,1);
                timeRhoHistogram    = zeros(numSim,1);
                timeRhoScott        = zeros(numSim,1);
                t                   = 1;
                
                save(fileName1);
            end
            
            angles = pi*(0:numAngles-1)/numAngles;
            angles = repmat(angles,1, ceil(numMeasurements(k)/numAngles));
            angles = angles(1:numMeasurements(k))';
            
           while t <= numSim,
                fprintf(['Iteraction ',num2str(t), '\n']);
                
                tic;
                S   = init_tables(Mph(i));
                psi = generate_cat_vector(alpha, phase, S);
                Rho = apply_loss(psi,etaState,S);
                
                Samples = homodyne_samples(-7,7,etaDetector,angles,Rho,S);
                
                % Constructs RhoML2 using the combined_optimization
                tic;
                [RhoML2, Diagnostics] = combined_optimization( Samples, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeML2(t) = toc;
                
                % Calculates the fidelity of RhoML2 and Rho
                fML2(t) = fidelity(RhoML2, Rho);
                    
                % Constructs the MHistogram matrix using the BinWidth (case option= 7 in the matrix_histogram function)
                tic;
                MHistogram = matrix_histogram2(Samples, 7, deltaq(j));  
                timeMhistogram(t) = toc;
                
                % Constructs the MScott array using the optimal width method (case option = 2 in the matrix_histogram function)
                tic;
                MScott = matrix_histogram(Samples, 2);
                timeMScott(t) =toc;
                
                [RhoHistogram, Diagnostics] = combined_optimization( MHistogram, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeRhoHistogram(t) = toc;
                
                % Constructs RhoScott using the combined_optimization
                tic;
                [RhoScott, Diagnostics] = combined_optimization( MScott, S, etaDetector, 0, maxIterations, stoppingCriterion);
                timeRhoScott(t) = toc;
                
                fHistogram(t)       = fidelity(RhoHistogram, Rho); % Calculates the fidelity of RhoHistogram and Rho
                fScott(t)           = fidelity(RhoScott, Rho);     % Calculates the fidelity of RhoScott and Rho
                fHistogramPsi(t)    = fidelity(RhoHistogram, psi); % Calculates the fidelity of RhoHistogram and Psi
                fML2Psi(t)          = fidelity(RhoML2, psi);       % Calculates the fidelity of RhoML2 and Psi
                fScottPsi(t)        = fidelity(RhoScott, psi);     % Calculates the fidelity of RhoScott and Rho
                
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
                meanFScottPsi           = mean(fScottPsi);     
                meanTimeRhoML2          = mean(timeML2);
                meanTimeMHistogram      = mean(timeMhistogram);
                meanTimeMScott          = mean(timeMScott);
                meanTimeRhoHistogram    = mean(timeRhoHistogram);
                meanTimeRhoScott        = mean(timeRhoScott);
                stdFHistogram           = std(fHistogram);        % Standard deviation of fHistogram
                stdFML2                 = std(fML2);              % Standard deviation of fML2
                stdFScott               = std(fScott);            % Standard deviation of fScott
                stdTML2                 = std(timeML2);
                stdTMHistogram          = std(timeMhistogram);
                stdTMScott              = std(timeMScott);
                stdTRhoHistogram        = std(timeRhoHistogram);
                stdTRhoScott            = std(timeRhoScott);
                
                home;
                fprintf('>> Progress: %.2f%%\n', t/numSim*100);
                t=t+1;
                
                save(fileName1);
            end
        end
    end
end