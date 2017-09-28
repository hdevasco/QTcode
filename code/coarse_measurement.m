function [H_int] = coarse_measurement(paramList,eta, S, returnMatrix)
% coarse_measurement returns the integrated measurement matrix along the bin where the integration
% limits x = a and x = b are the edges of the bin.

theta = paramList(1);
a = paramList(2);
b = paramList(3);
H_int = integral(@(x) homodyne_loss_measurement([x,theta],eta,S,returnMatrix),a,b,'ArrayValued',true);

