function [rho, Diagnostics] = rrhor_optimization(Measurements, photons, eta, linearForm, maxIterations, stop, rho, timeStart, continueOpt)
% performs r*rho*r optimization algorithm and return density matrix.
%
%   rho = rrhor_optimization(measurements, photons, eta, linearForm,
%   maxIterations, stop) estimates state rho using the homodyne measurement
%   data in Measurements.  Measurements can be a N-by-2 or N-by-3 array
%   containing N measurement results, where measurementArray(n,:) = [phase
%   angle, quadrature observed, number of observations], but if the third
%   column is absent it assumes each result was observed once.  Also,
%   Measurements can be the structre containing POVMs made by
%   make_measurement_struct.  photons is the number of photons that will be
%   included in the estimate, or it may be the structure made by
%   init_tables. eta is the efficiency of the homodyne detector. linearForm
%   may be a matrix whose dimension is equal to the dimension of the
%   density matrix rho.  rrhor_optimization will maximize L(rho) +
%   trace(linearForm*rho), where L is the loglikelihood.  If you only want
%   want to maximize likelihood, you can enter scalar 0 for the linearForm.
%   Each iteration will apply r*rho*r until maxIterations is reached or the
%   stopping criterion < stop.  The argument rho can be the density matrix
%   to initialize the iterations.  If rho is [], then the mixed state will
%   be used.  timeStart is an optional input that allows you to specify a
%   time at which iterations began. continue_opt is an optional input,
%   which if present prevents printing a warning when the iterations have
%   not converged.
%   
%   [rho, Diagnostics] = ... also returns a structure containing diagnostic
%   information about the progress of the maximization

if isnumeric(photons)
    S = init_tables(photons);
elseif isstruct(photons)
    S = photons;
end

Diagnostics.rhoArray = zeros(S.dimHilbertSpace,S.dimHilbertSpace,maxIterations+1);
Diagnostics.loglikelihoodList = zeros(maxIterations+1,1);
Diagnostics.objectiveList = zeros(maxIterations+1,1);
Diagnostics.linearTermList = zeros(maxIterations+1,1);
Diagnostics.rArrayList = zeros(S.dimHilbertSpace, S.dimHilbertSpace, maxIterations+1);
Diagnostics.stopList = zeros(maxIterations+1, 1);
Diagnostics.timeList = zeros(maxIterations+1, 1);
Diagnostics.ratio = zeros(maxIterations+1,1);

if isnumeric(Measurements)
    Measurements = make_measurement_struct(Measurements, eta, S);
end


% initial guess to begin maxlik r*rho*r algorithm, the maximally mixed
% state
if ~exist('rho', 'var') || (exist('rho', 'var') && isempty(rho)) || strcmp(rho, 'maxMix')
    rho = eye(S.dimHilbertSpace) ./ (S.dimHilbertSpace);
end

if not(exist('timeStart', 'var'))
    timeStart = tic;
end

% do iterations of r*rho*r
nIteration = 0;
R = make_r_struct(rho, Measurements, linearForm);
Diagnostics = add_diagnostics(R, Diagnostics, nIteration, timeStart);
while R.stop  > stop && nIteration < maxIterations
   nIteration = nIteration+1;
   rho = R.r * rho * R.r;
   rho= rho ./ trace(rho);
   R = make_r_struct(rho, Measurements, linearForm);
   Diagnostics = add_diagnostics(R, Diagnostics, nIteration, timeStart);
end
if nIteration == maxIterations && exist('continueOpt', 'var')~=1;
    warning('Tomography:fewIterations','maximum number of iterations reached before fidelity converged')
end

% trim extra rows from Diangostics
extraRow = nIteration + 2;
Diagnostics.rhoArray(:,:,extraRow:end) = [];
Diagnostics.loglikelihoodList(extraRow:end) = [];
Diagnostics.objectiveList(extraRow:end) = [];
Diagnostics.linearTermList(extraRow:end) = [];
Diagnostics.rArrayList(:,:,extraRow:end) = [];
Diagnostics.stopList(extraRow:end) = [];
Diagnostics.Measurements = Measurements;
Diagnostics.timeList(extraRow:end) = [];    
Diagnostics.typeList = ['start'; cellstr(repmat('RrhoR',nIteration, 1))];

Diagnostics.ratio (extraRow:end)=[];
rho = R.rho;

end

function Diagnostics = add_diagnostics(R, Diagnostics, nIteration, timeStart)
nRow = nIteration + 1;
Diagnostics.rhoArray(:,:,nRow) = R.rho;
Diagnostics.rArrayList(:,:,nRow) = R.r;
Diagnostics.stopList(nRow,1) = R.stop;
Diagnostics.loglikelihoodList(nRow,1) = R.loglike;
Diagnostics.objectiveList(nRow, 1) = R.objective;
Diagnostics.linearTermList(nRow, 1) = R.linearTerm;
Diagnostics.timeList(nRow,1) = toc(timeStart);
%Diagnostics.timeList(nRow,1) =0;

end