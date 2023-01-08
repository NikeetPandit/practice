function kernel = kernel_construct_notch(D0, Uk, Vk, M, N, type,varargin)
%%%% Inputs
% DO is cut off frequency 
% Uk, Vk are notch pairs (type array)
% M and N are dimensions of kernel size
% type is class string Options are: Gaussian, butterworth
% Vargin(1) specified butterworth order (type double)

%%%% Function Sumary
%D is is the distance between point (u,v) in frequency domian to notch of
%PxQ frequency rectangle (Digital Image Processing 4ed)
%
%Written by Nikeet Pandit

if length(varargin) ==1
    n = varargin{1}; 
end

if length(Uk) ~= length(Vk)
    disp('Notch pair Uk and Vk must be same length')
    return
end

type = lower(type);
u = 0:M-1; %frequency components
v = 0:N-1; 
[U, V] = meshgrid(u,v); 

if type == "gaussian" 
    kernel_mult = zeros(M,N,length(Uk)); 
    
    for i = 1:length(Uk) %if more than 1 notch frequencies specified
        
        % calculating positive notch distance mat. and positive notch kernel
        D =  my_calc_dist(U,V,Uk(i),Vk(i),M,N); 
        kernel_pos = fftshift(1-(exp(-(D.^2)./(2*(D0^2)))))'; %high pass kernel centered around notch Uk, Vk
        
        % calculating negative notch distance and negative notch kernel
        D =  my_calc_dist(U,V,-Uk(i),-Vk(i),M,N); 
        kernel_neg = fftshift(1-(exp(-(D.^2)./(2*(D0^2)))))'; %high pass kernel centered around notch -Uk, -Vk
        
        % notch filter are constructed as prod. of high pass kernels (for positive and negative notch)
        kernel_mult(:,:,i) = kernel_pos.*kernel_neg;  %symmetric kernel (about origin) (uncentered)
    end
    
    %Elementwise matrix along 3rd dimension to construct kernel...
    %if more than 1 positive and negative notch specified by user
    kernel = prod(kernel_mult,3); 
    
elseif type == "butterworth"
    kernel_mult = zeros(M,N,length(Uk)); 
    
    for i = 1:length(Uk)
        %Same as above except with butterworth kernel
        D =  my_calc_dist(U,V,Uk(i),Vk(i),M,N); 
        NumExp = D; 
        DenExp = D0; 
        kernel_pos = fftshift(1./(1+((NumExp./DenExp).^(2*n))))'; 
        D =  my_calc_dist(U,V,-Uk(i),-Vk(i),M,N); 
        NumExp = D; 
        DenExp = D0; 
        kernel_neg = fftshift(1./(1+((NumExp./DenExp).^(2*n))))'; 
        kernel_mult(:,:,i) = kernel_pos.*kernel_neg; 
    end
    
    kernel = prod(kernel_mult,3); 
end
end

%--- Calculating distance matrix for trasnfer functions 
function D =  my_calc_dist(U,V,Uk,Vk,M,N)
    D = hypot(U - M/2 - Uk, V - N/2 - Vk);
end


