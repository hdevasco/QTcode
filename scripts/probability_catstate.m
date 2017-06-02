tic;
clc;
clear;


alpha = 1;
    
number_photons = 1:1:15;

 % Probability
for i= 1:15; 

    p = 0:1:15; 
    
 % Normalization
 N = (2+2*exp(-2*alpha^2));
 
 P(i)= ((2/(N*factorial(i)))*(exp(-alpha^2)*(alpha.^2)^(i)+exp(-alpha^2)*(-alpha^2)^(i)));
 
 
end

  plot(number_photons,P);       
                