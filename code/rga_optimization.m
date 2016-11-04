function [rho, Diagnostics] = rga_optimization(Measurements, photons, eta, linearForm, maxRgaIterations, stop, display, Diagnostics, timeStart)
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

lamratio = 1.1;

if not(exist('display','var'))
    display=true;
end
if not(exist('timeStart','var'))
    timeStart=tic;
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
    rho = Diagnostics.rhoArray(:,:,end);
    maxIterations = nIteration + maxRgaIterations;
    R = make_r_struct(rho, Measurements, linearForm);
    Diagnostics.rhoArray(:, :, maxIterations+1) = zeros(S.dimHilbertSpace, S.dimHilbertSpace);
    Diagnostics.loglikelihoodList(maxIterations+1, 1) = 0;
    Diagnostics.objectiveList(maxIterations+1, 1) = 0;
    Diagnostics.linearTermList(maxIterations+1, 1) = 0;
    Diagnostics.rArray(:, :, maxIterations+1) = zeros(S.dimHilbertSpace, S.dimHilbertSpace);
    Diagnostics.stopList(maxIterations+1, 1) = 0;
    Diagnostics.timeList(maxIterations, 1) = 0;   
else
    if exist('Diagnostics','var') && isnumeric(Diagnostics)
        rho = Diagnostics;
        clear('Diagnostics')
    else
        rho = eye(S.dimHilbertSpace) ./ S.dimHilbertSpace;
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
    R = make_r_struct(rho, Measurements, linearForm);
    if isfield(R,'go')
        go = R.go;
    else    
        newCounterRow = [nIteration, 0, 0, R.objective, R.stop];    
        Diagnostics = add_diagnostics(R, Diagnostics, newCounterRow, timeStart);
    end

end


format short g

if display
    disp('     iteration,   maxStepSize, stepSize,    objective,       stop')
end

stepSize = 1;
newCounterRow = [nIteration, 0, 0, R.objective, R.stop];
iterationCounter = newCounterRow;

if display
    disp(newCounterRow)
end

go = true;

while R.stop > stop && nIteration < maxIterations && go
    nIteration = nIteration + 1;
    reduceStepSize = true;
    RI = make_ri_struct(R, Measurements);
    lambdaStart = lamratio*max([RI.mMaxEigenvalue; 0]);      % used to be just max([RI.mMaxEigenvalue; 0])
    % This value is chosen for lambdaStart to ensure that (lambdaStartRatio*lambda)*I-M is
    % always invertable.
    while reduceStepSize && go
        lambdaR = lambdaStart;
        increasedLambdaR = true;
        while increasedLambdaR
         aRI = ari_of_lambda(lambdaR, RI);
         checkStepSize = aRI.'*aRI;
         newCounterRow = [nIteration, stepSize, checkStepSize,  R.objective, R.stop];
         if display
            disp(newCounterRow)
         end
            iterationCounter = [iterationCounter; newCounterRow];
            if checkStepSize <= stepSize
                increasedLambdaR = false;
            elseif checkStepSize > stepSize
                lambdaR = 2*lambdaR;               % 2*lambdaR
                
                if lambdaR == inf                          % stop if lamba is too big (no convergence!)
                    go = false;
                    increasedLambdaR = false;
                end   
            end
        end %finding lambdaR
        
        if go
            a = unvectorize_r_i(aRI);
            rhoTest = (RI.rhoSqrt+a)*(RI.rhoSqrt+a');
            rhoTest = rhoTest/trace(rhoTest);
            RTest = make_r_struct(rhoTest, Measurements, linearForm);
            if RTest.objective > R.objective
                reduceStepSize = false;
                R = RTest;
                newCounterRow = [nIteration, 0, 0, R.objective, R.stop];
                if display
                    disp(newCounterRow)
                end
                iterationCounter = [iterationCounter; newCounterRow];
            else
                stepSize = stepSize/2;
                reduceStepSize = true;
                if stepSize == 0
                    go = false;
                end
            end
        end

    end % end finding small enough stepSize
    Diagnostics = add_diagnostics(R, Diagnostics, nIteration, timeStart);
end % rga iterations

if display
    if nIteration == maxIterations
        warning('Tomography:fewIterations','maximum number of iterations reached before stopping criterion converged')
    end
end

if ~go 
    warning('Tomography:noConvergence','maximum step size converged to zero before stopping criterion converged!')
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

end

function Diagnostics = add_diagnostics(R, Diagnostics, counterRow, timeStart)

nIteration = counterRow(1);

nRow = nIteration + 1;
Diagnostics.rhoArray(:,:,nRow) = R.rho;
Diagnostics.rArray(:,:,nRow) = R.r;
Diagnostics.stopList(nRow,1) = R.stop;
Diagnostics.loglikelihoodList(nRow,1) = R.loglike;
Diagnostics.objectiveList(nRow, 1) = R.objective;
Diagnostics.linearTermList(nRow, 1) = R.linearTerm;
Diagnostics.timeList(nRow, 1) = toc(timeStart);
if nIteration == 0
    Diagnostics.typeList(nRow, 1) = cellstr('start');
else
    Diagnostics.typeList(nRow, 1) = cellstr('rga');
end


end

