function [rho, Diagnostics, nIteration] = rga_optimization_upgrade(Measurements, photons, eta, linearForm, maxRgaIterations, stop, display, Diagnostics, timeStart)

% performs iterations of the regularized gradient ascent optimization
%   rga_optimization(Measurements, photons, eta, linearForm, maxIterations,
%   stop, display, Diagnostics) performs iteraions of the regularized gradient ascent
%   maximization. It will maximize L(rho)+trace(linearForm*rho), where L is
%   the loglikelihood. measurements can be a N-by-2 or N-by-3 array
%   containing N measurement results, where measurementArray(n,:) = [phase
%   angle, quadrature observed, number of observations], but if the third
%   column is absent it assums each result was observed once.  Also,
%   measurements can be the structre containing POVMs made by
%   make_measurement_struct.  photons is the number of photons to include
%   in the reconstructed density matrix, or it may be the structure made by
%   init_tables.  eta is the efficiency of the homodyne detector. Each
%   iteration will apply the regularized gradient ascent until
%   maxRgaIterations is reached or the the stopping criterion < stop.
%   Diagnostics is optional and if present should contain the diagnostic
%   data created by previously running tomography, for example, starting
%   with rrhor_optimization and then switching to rga_optimization.  Also
%   Diagnostics may just be a density matrix, in which case it will be used
%   as the matrix to begin the iterations.  (Otherwise the mixed state is
%   used).  rga_iterations will continue where the previous effort left
%   off.
%
%   [rho, Diagnostics] = ...also returns a structure with diagnostic
%   information about the progress of the maximization.
%
%
%	This is the "upgraded" version of RGA, using the numerical optimization tricks 
%	from the textbook involving changing the size of the trust region at every step.  
%	The maximum step size can increase, decrease, or stay the same at each step depending on
% 	if the ratio between the predicted increase and the actual increase is large, very small, or in between.
%
%	Also, lambda starts at 1.1*maxEigvalue instead of the original 2

if not(exist('display','var'))
    display=true;
end

if isnumeric(photons)
    S = init_tables(photons);
elseif isstruct(photons)
    S = photons;
end

if isnumeric(Measurements)
    Measurements = make_measurement_struct(Measurements, eta, S);
end

if exist('Diagnostics', 'var') && isstruct(Diagnostics)
    nIteration = size(Diagnostics.loglikelihoodList, 1) - 1;
    maxIterations = nIteration + maxRgaIterations;
    rho = Diagnostics.rhoArray(:,:,end);
    R = make_r_struct(rho, Measurements, linearForm);
    Diagnostics.rhoArray(:, :, maxIterations+1) = zeros(S.dimHilbertSpace, S.dimHilbertSpace);
    Diagnostics.loglikelihoodList(maxIterations+1, 1) = 0;
    Diagnostics.objectiveList(maxIterations+1, 1) = 0;
    Diagnostics.linearTermList(maxIterations+1, 1) = 0;
    Diagnostics.rArray(:, :, maxIterations+1) = zeros(S.dimHilbertSpace, S.dimHilbertSpace);
    Diagnostics.stopList(maxIterations+1, 1) = 0;
    Diagnostics.timeList(maxIterations, 1) = 0;   
    
    Diagnostics.ratio(maxIterations, 1) =0;
else
    if exist('Diagnostics','var') && isnumeric(Diagnostics)
        rho = Diagnostics;
        clear('Diagnostics')
    else
        rho = eye(S.dimHilbertSpace) ./ S.dimHilbertSpace;
    end
    
    if not(exist('timeStart','var'))
        timeStart=tic;
    end    
    
    nIteration = 0;
    maxIterations = nIteration + maxRgaIterations;
    Diagnostics.rhoArray = zeros(S.dimHilbertSpace,S.dimHilbertSpace,maxIterations+1);
    Diagnostics.loglikelihoodList = zeros(maxIterations+1,1);
    Diagnostics.objectiveList = zeros(maxIterations+1,1);
    Diagnostics.linearTermList = zeros(maxIterations+1,1);
    Diagnostics.rArray = zeros(S.dimHilbertSpace, S.dimHilbertSpace, maxIterations+1);
    Diagnostics.stopList = zeros(maxIterations+1, 1);
    Diagnostics.timeList = zeros(maxIterations+1, 1);
    Diagnostics.ratio = zeros(maxIterations+1,1);
    
    R = make_r_struct(rho, Measurements, linearForm);
    
    newCounterRow = [nIteration, 0, 0, R.objective, R.stop, 0];
    
    Diagnostics = add_diagnostics(R, Diagnostics, newCounterRow, timeStart);
    
end




%------- ratio bounds for changing step size -
ratioUp=.75;        % upper bound for ratio - if ratio is greater than this, increse step size
ratioLo=.25;        % lower bound for ratio - if ratio is less than this, decrease step size
ratioAccept = .05;   % ratio must be at least this big in order to accept this step


%------- initialize step size-----------------
stepSize = 1;       % initial step size
maxStepSize =4;     % maximum step size

stepIncreaseRate = 2;       % increase step much by this rate if ratio is high enough
stepDecreaseRate = .5;      % decrease step much by this rate if ratio is low enough

lambdaStartRatio = 1.1;
lambdaIncreaseRate = 2; % sqrt(2)


%------- display counter ---------------------
format short g
increaseRatio=0;
if display
    disp('     iteration,   maxStepSize, stepSize,   objective,        stop,    increase ratio')
end
newCounterRow = [nIteration, 0, 0, R.objective, R.stop, increaseRatio];
if display   
    disp(newCounterRow)
end
iterationCounter = newCounterRow;

%---------------------------------------------------
%-- begin optimization -----------------------------
while R.stop > stop && nIteration < maxIterations
    nIteration = nIteration + 1;
    acceptStep = false;
    RI = make_ri_struct(R, Measurements);
    lambdaStart = lambdaStartRatio*max([RI.mMaxEigenvalue; 0]);
    % This value is chosen for lambdaStart to ensure that (lambdaStartRatio*lambda)*I-M is
    % always invertable.
    while not(acceptStep)
        lambdaR = lambdaStart;
        increasedLambdaR = true;
        while increasedLambdaR
            aRI = ari_of_lambda(lambdaR, RI);
            checkStepSize = aRI.'*aRI;
            newCounterRow = [nIteration, stepSize, checkStepSize,  R.objective, R.stop, increaseRatio];
            if display
               disp(newCounterRow)
            end
            iterationCounter = [iterationCounter; newCounterRow];
            if checkStepSize <= stepSize
                increasedLambdaR = false;
            elseif checkStepSize > stepSize
                lambdaR = lambdaIncreaseRate*lambdaR;
            end
        end %finding lambdaR

        a = unvectorize_r_i(aRI);
        rhoTest = (RI.rhoSqrt+a)*(RI.rhoSqrt+a');
        rhoTest = rhoTest/trace(rhoTest);
        RTest = make_r_struct(rhoTest, Measurements, linearForm);
        
        %------- calculate ratio of actual to projected increase
        
        actualIncrease = (RTest.objective)-(R.objective);
        predictIncrease = RI.v.'*aRI + 1/2*aRI.'*RI.m*aRI + trace_of_product(linearForm,(RTest.rho-R.rho));
        
        increaseRatio = actualIncrease/predictIncrease;
        
        %------- check ratio bounds to see if stepSize should be changed
        
        if increaseRatio > ratioUp
            stepSize = min(maxStepSize, stepIncreaseRate * stepSize);
        elseif increaseRatio < ratioLo
            stepSize = stepDecreaseRate * stepSize;
        end
        
        %------ check to see if step is accepted        

        if increaseRatio > ratioAccept
            acceptStep=true;                % step accepted
                        
            R = RTest;
            newCounterRow = [nIteration, 0, 0, R.objective, R.stop, increaseRatio];
            if display
                disp(newCounterRow) 
            end
            iterationCounter = [iterationCounter; newCounterRow];
         
        else    
            acceptStep=false;               % step not accepted
        end
        
    end % end finding small enough stepSize
    Diagnostics = add_diagnostics(R, Diagnostics, newCounterRow, timeStart);
end % rga iterations

if nIteration == maxIterations
    warning('Tomography:fewIterations','maximum number of iterations reached before fidelity converged')
end

% trim extra rows from Diangostics
extraRow = nIteration + 2;
Diagnostics.rhoArray(:,:,extraRow:end) = [];
Diagnostics.loglikelihoodList(extraRow:end) = [];
Diagnostics.linearTermList(extraRow:end) = [];
Diagnostics.objectiveList(extraRow:end) = [];
Diagnostics.rArray(:,:,extraRow:end) = [];
Diagnostics.stopList(extraRow:end) = [];
Diagnostics.Measurements = Measurements;
Diagnostics.iterationCounter = iterationCounter;
Diagnostics.typeList(:,:,extraRow:end) = [];
Diagnostics.timeList(extraRow:end) = [];
rho = R.rho;

Diagnostics.ratio(extraRow:end)=[];
end

function Diagnostics = add_diagnostics(R, Diagnostics, counterRow, timeStart)

nIteration = counterRow(1);
increaseRatio = counterRow(6);

nRow = nIteration + 1;
Diagnostics.rhoArray(:,:,nRow) = R.rho;
Diagnostics.rArray(:,:,nRow) = R.r;
Diagnostics.stopList(nRow,1) = R.stop;
Diagnostics.loglikelihoodList(nRow,1) = R.loglike;
Diagnostics.objectiveList(nRow, 1) = R.objective;
Diagnostics.linearTermList(nRow, 1) = R.linearTerm;
Diagnostics.timeList(nRow, 1) = toc(timeStart);
Diagnostics.increaseRatio(nRow,1)=increaseRatio;

if nIteration == 0
    Diagnostics.typeList(nRow, 1) = cellstr('start');
else
    Diagnostics.typeList(nRow, 1) = cellstr('rga');
end    

end


