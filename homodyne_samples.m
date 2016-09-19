 function [samples] = homodyne_samples(min_q, max_q, eta, angles, rho, S)
%HOMODYNE_SAMPLES simulates homodyne measurement of state rho.
%
%   SAMPLES=HOMODYNE_SAMPLES(MIN_Q, MAX_Q, ETA, ANGLES, RHO, S) generates
%   simulated homodyne measurements of the input state RHO.  It only
%   consideres possible measurement quadrature results between MIN_Q and
%   MAX_Q.  ETA is the efficiency of the homodyne detection.  ANGLES is a
%   vector that contains a list of phase angles at which to make the
%   measurements.  If you want to measure a particular angle several times,
%   that angle must be listed several times in the vector ANGLES.  S the
%   structre created by init_tables.  HOMODYNE_SAMPLES returns an N-by-3
%   array.  N = length(ANGLES).  Each row contains [phase angle, quadrature
%   result, 1].  The 1 means that result was observed 1 time.  It is
%   included for compatibilty with repeated measurement results.

number_of_angles = size(angles,1);
samples = [angles,zeros(number_of_angles,1),ones(number_of_angles,1)];

n=S.photons;

if eta ~= 1
    rho = apply_loss(rho,eta,S);
end

% This max_pdf is the maximum of the pdf for a pure squeezed state with mean
% number of photons equal to mean_n.  I don't know why this wasn't working,
% so I multiplied max_pdf by 2.
mean_n = mean_photons(rho);
max_pdf = exp(asinh(mean_n.^2))/sqrt(pi)*2;

% Choose random measurement results q according to probability distribution
% given by trace(rho*POVM(q,theta)).
for j = 1:number_of_angles
    pdfguess = 1; % a value which is greater than pdf_at_qguess in order to begin the loop
    pdf_at_qguess = 0;
    while(pdf_at_qguess < pdfguess)
        qguess = rand * (max_q - min_q) + min_q;
        quadrature_vector_at_qguess = homodyne_loss_measurement([angles(j);qguess], 1, S,'return_vector');
        pdfguess = rand * max_pdf;
        if isvector(rho)
            pdf_at_qguess = abs(quadrature_vector_at_qguess'*rho)^2;
        else
            pdf_at_qguess = quadrature_vector_at_qguess'*rho*quadrature_vector_at_qguess;
        end
        if pdf_at_qguess > max_pdf
            error('Tomography:max_pdfTooSmall','Maximum pdf value is too small.')
        end
    end
    samples(j,2) = qguess;
end