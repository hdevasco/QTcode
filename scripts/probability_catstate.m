tic;
clc;
clear;

% Probability distribution of n photons in cat state

alpha = 1;
    
number_photons = 0:1:16;

 % Probability
for i= 1:17; 

    
 % Normalization
 N = (2+2*exp(-2*alpha^2));
 
 P(i)= ((2/(N*factorial(i-1)))*(exp(-alpha^2)*(alpha^2)^(i-1)+exp(-alpha^2)*(-alpha^2)^(i-1)));
 
 
end

  plot(number_photons,P);       
                