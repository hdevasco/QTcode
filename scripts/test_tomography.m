tic;
clc;
clear;


%Maximum number of photons

maxPhotonNumber = [10];

% Number of measurements
numMeasurements    = 20000;

option = [100]

%deltaq = [0.3];


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
        for j=1:length(option),
            
            fileName1 = ['maxPhNum',num2str(maxPhotonNumber(i)),'option',num2str(option(j)),'alpha',num2str(alpha(k)),'.mat'];
            
            if exist(fileName1,'file') == 2,
                load(fileName1);
            else
                fML2                = zeros(numSim,1);
                fScott              = zeros(numSim,1);
                fHistogram          = zeros(numSim,1);
                fHistogram2          = zeros(numSim,1);
                fHistogramPsi       = zeros(numSim,1);
                fML2Psi             = zeros(numSim,1);
                fScottPsi           = zeros(numSim,1);
                fidelityDiff        = zeros(numSim,1);
                fidelityDiff2       = zeros(numSim,1);
                t                   = 1;
                
                save(fileName1);
            end
            
            while t <= numSim,
                fprintf(['Iteraction ',num2str(t), '\n']);
                
                angles = pi*(0:numAngles-1)/numAngles;
                angles = repmat(angles,1, ceil(numMeasurements/numAngles));
                angles = angles(1:numMeasurements)';
                
                %generate state
                
                S   = init_tables(maxPhotonNumber(i));
                psi = generate_cat_vector(alpha(k), phase, S);
                Rho = apply_loss(psi,etaState,S);
                
                
                Samples = homodyne_samples(-7,7,etaDetector,angles,Rho,S);
                
                
                [RhoML2, Diagnostics] = combined_optimization( Samples, S, etaDetector, 0, maxIterations, stoppingCriterion);
                
                fML2(t) = fidelity(RhoML2, Rho);
                
                MHistogram = matrix_histogram(numAngles, Samples, option(j),'integral');
                
                
                [RhoHistogram, Diagnostics] = combined_optimization( MHistogram, S, etaDetector, 0, maxIterations, stoppingCriterion);
                
                MHistogram2 = matrix_histogram(numAngles, Samples, option(j),'integral');
                
                
                [RhoHistogram2, Diagnostics] = combined_optimization( MHistogram2, S, etaDetector, 0, maxIterations, stoppingCriterion);
                
                MScott = matrix_histogram(numAngles, Samples, 8,'integral');
                
                
                [RhoScott, Diagnostics] = combined_optimization( MScott, S, etaDetector, 0, maxIterations, stoppingCriterion);
                
                MScott2 = matrix_histogram_new(numAngles, Samples, 8,'integral');
                
                
                [RhoScott2, Diagnostics] = combined_optimization( MScott2, S, etaDetector, 0, maxIterations, stoppingCriterion);
                
                fHistogram(t)       = fidelity(RhoHistogram, Rho);
                fHistogram2(t)       = fidelity(RhoHistogram2, Rho);
                fScott(t)           = fidelity(RhoScott, Rho);     % Calculates the fidelity of RhoScott and Rho
                fScott2(t)          = fidelity(RhoScott2, Rho);
                fidelityDiff(t)     = fML2(t)-fHistogram(t);       % Difference between fidelities [fidelity(RhoML2,Rho)-fidelity(RhoHistogram,Rho)]
                fidelityDiff2(t)    = fML2(t)-fScott(t);           % Difference between fidelities [fidelity(RhoML2,Rho)-fidelity(RhoScott,Rho)]
                wHistogram          = mean(fidelityDiff);
                wHistogram2         = mean(fidelityDiff2);
                Dhistogram          = std(fidelityDiff);           % Standard deviation of the difference between the fidelities(fML2,fHistogram)
                Dhistogram2          = std(fidelityDiff2);         % Standard deviation of the difference between the fidelities(fML2,fScott)
                meanFML2            = mean(fML2);
                meanFHistogram      = mean(fHistogram);
                meanFHistogram2     = mean(fHistogram2);
                meanFScott          = mean(fScott);
                meanFScott2         = mean(fScott2);
                meanFML2Psi         = mean(fML2Psi);
                meanFHistogramPsi   = mean(fHistogramPsi);
                meanFScottPsi           = mean(fScottPsi);
                stdFHistogram           = std(fHistogram);        % Standard deviation of fHistogram
                stdFHistogram2           = std(fHistogram2);
                stdFML2                 = std(fML2);              % Standard deviation of fML2
                stdFScott               = std(fScott);
                stdFScott2               = std(fScott2);            % Standard deviation of fScott
                home;
                fprintf('>> Progress: %.2f%%\n', t/numSim*100);
                t=t+1;
                
                save(fileName1);
            end
        end
    end
end