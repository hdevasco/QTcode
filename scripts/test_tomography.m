clc;
clear;

savefile = 'raw_quadratures.mat';
% Number of angles equally spaced from 0 to pi
numAngles = 20;

if exist(savefile, 'file')~=2
    
    % Maximum number of photons
    maxPhotonNumber = 10;

    % Number of measurements
    numMeasurements    = 20000;

    etaDetector        = 0.9;
    maxIterations      = 2000;
    stoppingCriterion  = 0.01;
    alpha              = 1;
    phase              = 0;
    etaState           = 0.95;

    angles = pi*(0:numAngles-1)/numAngles;
    angles = repmat(angles,1, ceil(numMeasurements/numAngles));
    angles = angles(1:numMeasurements)';

    % Generate state
    S   = init_tables(maxPhotonNumber);
    psi = generate_cat_vector(alpha, phase, S);
    Rho = apply_loss(psi,etaState,S);

    Samples = homodyne_samples(-7,7,etaDetector,angles,Rho,S);

    save(savefile, 'Samples')

elseif exist(savefile,'file')==2
    load(savefile)
end


deltaq    = 0.5;

optList   = [1, 2, 3, 4, 5, 6, 7, 8, 100];
tic;
for j = 1:length(optList)
    MHistogram_old_center = matrix_histogram(numAngles, Samples, optList(j),'center', deltaq);
    MHistogram_new_center = matrix_histogram_new(numAngles, Samples, optList(j),'center', deltaq);
    if ~isequal(MHistogram_old_center, MHistogram_new_center)
        error('Error when option is %d and center.', optList(j))
    end
    
    MHistogram_old_int = matrix_histogram(numAngles, Samples, optList(j),'integral', deltaq);
    MHistogram_new_int = matrix_histogram_new(numAngles, Samples, optList(j),'integral', deltaq);
    if ~isequal(MHistogram_old_int, MHistogram_new_int)
        error('Error when option is %d and integral.', optList(j))
    end
    
end
toc
sprintf('Tudo bem!')
