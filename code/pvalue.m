function [pStruct] = pvalue(measurementData1, measurementData2, photons, linearForm, eta, pTarget, maxIteration, stop, numTimesSimulate)

%
% Split measurement data into 2 halves, measurementData1 and
% measurementData2.  Reconstruct rho_Ml from the first data set.
% Then numTimesSimulate*nTotalMeasurements measurements
% are simulated using rho_Ml.  This third data set is used to find an
% approximation for the density matrix sigma that minimizes the
% KL-divergence to rho_Ml.
%
%       photons - should be either a number or a structure from init_tables
%       linearForm - the linear functional 
%       eta - efficiency for the input data sets
%       pTarget - the target significance level (I used 0.05)\
%       maxIteration, stop - parameters for the reconstructions done inside
%               of this routine
%       numTimesSimulate - amount of data to simulate for KL-divergence
%               minimization
%
%

maxIteration1 = maxIteration;
maxIteration2 = maxIteration;
stop1 = stop;
stop2 = stop;


if isstruct(photons)
    s = photons;
else
    s = init_tables(photons);
end

Data1Struct = make_measurement_struct(measurementData1, eta, s );
[rhoML, ~] = combined_optimization(Data1Struct, s, eta, 0, maxIteration1, stop1, 'maxMix', 0, 0);

numDataSimulated = numTimesSimulate*Data1Struct.nTotalMeasurements;
angles = rand(size(1:numDataSimulated))'*2*pi;
[DataML] = homodyne_samples(-12, 12, 1, angles, rhoML, s);
DataMLStruct = make_measurement_struct(DataML, 1, s );

Data2Struct = make_measurement_struct(measurementData2, eta, s );

increaseLambda = true;
lambda = 2^(-2);
%lambda = 1;
counter = 0;

go = true;

likRhoMlD2 = 1;   % likelihood of rhoML (from data set 1) given data set 2
for n = 1:Data2Struct.nMeasurements
    trRhoMlD2 = trace_of_product(rhoML, Data2Struct.povmArray(:,:,n)) ;
    likRhoMlD2 = likRhoMlD2*((trRhoMlD2)^ Data2Struct.measurementArray(n,3));
end

while increaseLambda && go
    counter = counter + 1;
    lambdas(counter) = lambda;
    
    increaseLambda = false;
    if lambdas(counter) < 10000
            [sigmaLambda, ~] =  combined_optimization(DataMLStruct, s, eta, lambda*linearForm, maxIteration2, stop2, rhoML, 0, 0);
    else
        go = false;
    end
    
    statistic = likRhoMlD2;
    for n = 1:Data2Struct.nMeasurements
        trSigmaD2 = trace_of_product(sigmaLambda, Data2Struct.povmArray(:,:,n)) ;
        statistic = statistic*((1/trSigmaD2)^ Data2Struct.measurementArray(n,3));
    end
    pLambda = real(1/statistic);
    f = trace_of_product(sigmaLambda, linearForm);
    fs(counter)=f;
    ps(counter) = pLambda;
    if pLambda > pTarget
        increaseLambda = true;
        lambda = 1.1892*lambda;  % 2^(1/4)
    else
        %lambda
    end
    
end

f = trace_of_product(sigmaLambda, linearForm);
pValue = pLambda;


pStruct = struct('pvalueList', ps, 'fList', fs, 'rhoML', rhoML, 'lambdas', lambdas);



