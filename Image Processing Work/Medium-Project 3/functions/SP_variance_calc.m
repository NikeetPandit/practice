%--- Function to calculate S&P noise variance
% 
%Written by Nikeet Pandit

function SPvar = SP_variance_calc(noise_density)
Ps = noise_density/2; Pp = Ps; 
K = 2; %intensity levels 
SPmean = 0*Pp + K*(1-Ps-Pp) + (2^K-1)*Ps; %from Digital Image Processing 4ed
SPvar = (0-SPmean)^2*Pp + (K-SPmean)^2*(1-Ps-Pp) + (2^K-1)^2*Ps; 
end

