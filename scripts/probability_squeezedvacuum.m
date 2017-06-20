tic;
clc;
clear;

% Probability distribution of n photons in squeezed state

r = 0.14;

n = 10
    
alpha = 0:1:8;

theta = 0;

for i= 1:9; 

 z = ((i-1)+(i-1)*exp(1j*theta)*tanh(r))/((2*exp(1j*theta)*tanh(r))^(-1/2));
 
 
 t = (-((i-1)^2)-(1/2)*(((i-1)^2)*exp(1j*theta)+((i-1)^2)*exp(-1j*theta))*tanh(r));
 
 
 P(i)= (((1/2)*(tanh(r))^n)/((factorial(n)*cosh(r)))*exp(t)*abs(hermiteH(10,z))^2);
 
 
end

  plot(alpha,P);       
                