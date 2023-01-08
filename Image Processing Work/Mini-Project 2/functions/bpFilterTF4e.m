function kernel = bpFilterTF4e(type, M, N, C0, W, n)
%%%% Inputs
% CO is centre of band
% W is width of the band
% M and N are dimensions of kernel size
% type is class string 
% Options are: Gaussian, butterworth, ideal

%%%% Function Sumary
% Calculates centered frequency domain filtering kernel
%
% Written by Nikeet Pandit

kernel = 1 - brFilterTF4e(type, M, N, C0, W, n); %Bandpass filter = 1 - Bandreject filter 
end
