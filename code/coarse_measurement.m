function [H_int] = coarse_measurement(paramList,edges)
% coarse_measurement returns the integrated measurement matrix along the bin where the integration
% limits x = a and x = b are the edges of the bin.

for i=1:length(edges)-1,
   
H_int = integral(@(x) homodyne_loss_measurement([x,2],eta,S,'return matrix'),edges(i),edges(i+1),'ArrayValued',true);

end

