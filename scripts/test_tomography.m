tic;
clc;
clear;
%Maximum number of photons

maxPhotonNumber = [10];

% Number of measurements
numMeasurements    = 20000;

option = [100];

%deltaq = [0.5];

% Number of simulations
numSim             = 1;

etaDetector        = 0.9;
maxIterations      = 2000;
stoppingCriterion  = 0.01;
alpha              = [1] ;
phase              = 0;
etaState           = 0.95;

% Number of angles equally spaced from 0 to pi
numAngles          = 20;


for i=1:length(maxPhotonNumber),
    for k=1:length(alpha),
        for j=1:length(option),
            
            fileName1 = ['maxPhNum',num2str(maxPhotonNumber(i)),'option',num2str(option(j)),'alpha',num2str(alpha(k)),'.mat'];
            
            if exist(fileName1,'file') == 2,
                load(fileName1);
            else
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
                
                MHistogram_old = matrix_histogram(numAngles, Samples, option(j),'center');
                
                MHistogram_new = matrix_histogram_new(numAngles, Samples, option(j),'center');
                
                isequal(MHistogram_old, MHistogram_new)
       
                home;
                fprintf('>> Progress: %.2f%%\n', t/numSim*100);
                t=t+1;
                
                save(fileName1);
            end
        end
    end
end