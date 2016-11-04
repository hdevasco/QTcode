function [ rho, Diagnostics ] = combined_optimization( Measurements, photons, eta, linearForm, maxRgaIterations, stop, rho, display, upgrade)
% performs optimization using both RrhoR and gradient ascent algorithms
%   rho = combined_iterations(Measurements, photons, eta, linearForm,
%   maxRgaIterations, stop, rho) performs optimization over density
%   matrices using first R*rho*R iterations, then switches to the
%   regularized gradient ascent. It will maximize
%   L(rho)+trace(linearForm*rho).  Measurements can be a N-by-2 or N-by-3
%   array containing N measurement results, where measurementArray(n,:) =
%   [phase angle, quadrature observed, number of observations], but if the
%   third column is absent it assums each result was observed once.  Also,
%   Measurements can be the structre containing POVM elements made by
%   make_measurement_struct.  photons is the number of photons to include
%   in the reconstructed density matrix, or photons can be the structure
%   made by init_tables.  eta is the efficiency of the homodyne detector.
%   It begins by performing several iterations of R*rho*R, then performs
%   maxRgaIterations, but it will stop if stopping criterion > stop.  It
%   returns the maximum likelihood density matrix estimate, rho.
%
%   rho = combined_iterations(..., rho) uses rho as the initial density.
%   It defaults to the maximally mixed state, which can also be specified
%   with 'maxMix'.
%
%   display is a logical that toggles whether to display steps during
%   RGA optimization
%
%   upgrade is a logical that toggles whether to use rga_optimization
%   or rga_optimization_upgrade
%
%   [rho, Diagnostics] = combined_iterations(...) also returns a structure
%   containing diagnostic information about the algorithms' progress.

if ~exist('display','var')
    display = 0;
end
if ~exist('upgrade','var')
    upgrade = 1;
end
if isnumeric(photons)
    S = init_tables(photons);
elseif isstruct(photons)
    S = photons;
end
if isnumeric(Measurements)
    Measurements = make_measurement_struct(Measurements, eta, S);
end
if ~exist('rho', 'var') || isequal(rho, 'maxMix') 
    rho = eye(S.dimHilbertSpace)/S.dimHilbertSpace;
end

% TODO: there is a bug when maxIterations < initialIterations!!

initialIterations = round(S.dimHilbertSpace^2/4);
% We will do this number of R*rho*R iterations before switching to RGA.
% This takes usually about the same amount of time that one iteration of
% RGA would take.

timeStart=tic; % Before I added this line, time was not kept track of properly in combined_optimization. 
		% Now, timestamps can be used as inputs into the optimization algorithms -mgirard

[~, Diagnostics] = rrhor_optimization(Measurements, S, eta, linearForm, initialIterations, stop, rho, timeStart, 'continue0pt');

if upgrade == 1
    [rho, Diagnostics] = rga_optimization_upgrade(Measurements, S, eta, linearForm, maxRgaIterations, stop, display, Diagnostics, timeStart);
else
    [rho, Diagnostics] = rga_optimization(Measurements, S, eta, linearForm, maxRgaIterations, stop, display, Diagnostics, timeStart);
end

end

