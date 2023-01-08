%%% Function add noise
% Inputs 
% (1) Image (input)
% (2) [a, b] inputs for Gauss
% (3) S&P density
% (4) [a, b] for uniform noise

% Outputs
% (1) Original Image
% (2) Gauss noise added image
% (3) S&P added image
% (4) Uniform noise added image
%
%Written by Nikeet Pandit


function [Im, Im_SP, Im_Gaus, Im_Uform] = add_noise_func(Im, ab_gauss, d_sp, ab_uform)
[M, N] = size(Im);     
Im_Gaus = imnoise(Im,'gaussian',ab_gauss(1),ab_gauss(2));
Im_SP = imnoise(Im,'salt & pepper',d_sp); 
Im_Uform = Im+ ab_uform(1) + (ab_uform(2) - ab_uform(1))*rand(M,N); %uniform noise
end
