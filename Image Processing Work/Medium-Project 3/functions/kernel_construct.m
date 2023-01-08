function kernel = kernel_construct(D0, P, Q, type,varargin)
%%%% Inputs
% DO is cut off frequency 
% P and Q are dimensions of padded kernel size
% type is class string 
% Options are: Gaussian, ...

%%%% Function Sumary
%D is is the distance between point (u,v) in frequency domian at centre of
%PxQ frequency rectangle (Digital Image Processing 4ed)
%
%Written by Nikeet Pandit

if length(varargin) ==1
    n = varargin{1}; 
end

type = lower(type);
u = 0:P-1; %frequency components
v = 0:Q-1; 
[U, V] = meshgrid(u,v); 
D = hypot(U-P/2, V-Q/2); %calculating distances

if type == "gaussian"
    kernel = fftshift(exp(-(D.^2)./(2*(D0^2))))'; %Shifts kernel to nominal fft2 position (uncentered)
    
elseif type == "butterworth"
    NumExp = D; 
    DenExp = D0; 
    kernel = fftshift(1./(1+((NumExp./DenExp).^(2*n))))'; %Butterworth T.F. (uncentered)
  
end