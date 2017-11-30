function [H_int] = coarse_measurement(paramList,eta, S, returnMatrix)

% coarse_measurement returns the homodyne_loss_measurement.m measurement matrix integrated in 
% the integration interval x = a to x = b, where "a" and "b" are the edges of the histogram bin.

theta = paramList(1);
a = paramList(2);
b = paramList(3);
% matlab works with standard error parameters for integral function in:
% 'AbsTol' = 1e-10
% 'RelTol' = 1e-6
% Soon I increased these parameters to:
% 'AbsTol' = 1e-6
% 'RelTol' = 1e-2

 H_int = integral(@(x) homodyne_loss_measurement([theta,x],eta,S,returnMatrix),a,b,...
 'RelTol' , 1e-2, 'AbsTol' , 1e-6, 'ArrayValued',true);

 %H_int = integral(@(x) homodyne_loss_measurement([theta,x],eta,S,returnMatrix),a,b, 'ArrayValued',true);

end