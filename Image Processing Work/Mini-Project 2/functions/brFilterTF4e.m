function kernel = brFilterTF4e(type, M, N, C0, W, n)
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

type = lower(type);

%% Calculate D(u,v)
%D is is the dist. between point (u,v) in frequency at centre of MxN freq. rectangle

u = 0:M-1; %frequency components
v = 0:N-1; 
[U, V] = meshgrid(u,v); 
D = hypot(U-M/2, V-N/2); %calculating distances


if type == "gaussian" 
    NumExp = D.^2 - C0^2; 
    DenExp = D*W;
    kernel = 1 - exp(-((NumExp./DenExp).^2))'; %Gaussian T.F.
    
elseif type == "butterworth"
    NumExp = D*W; 
    DenExp = D.^2-C0^2; 
    kernel = (1./(1+((NumExp./DenExp).^(2*n))))'; %Butterworth T.F.
        
elseif type == "ideal"
    Logic_False_Cond = C0 - (W/2) <= D & D <= C0 + (W/2); %Ideal T.F.
    kernel = double(~Logic_False_Cond);
    
else
    error("Invalid Selection for type"); 
    
end