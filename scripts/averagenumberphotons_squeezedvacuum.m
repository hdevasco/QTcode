% (<|alpha,r|n|alpha,r>) Average number of protons according to the alpha value (Squeezed State)

average_number_photons = zeros(9,1);

r= 0.14;

alpha = [0,1, 2, 3, 4, 5, 6, 7, 8];

theta = 0;

for j=1:1:9;


average_number_photons(j)= (j-1)^2*(cosh(r)^2+sinh(r)^2)-(j-1)^2*(exp(1i*theta)*sinh(r)*cosh(r))-(j-1)^2*(exp(-1i*theta)*sinh(r)*cosh(r))+sinh(r)^2;

end


plot(alpha, average_number_photons)