%% Function low passes filters an image 
%In function DFT may mean DFT or DCT!!!

% Inputs: 
% (1) Im -image file read in by Imread(File)
% (2) D0 -cut off frequency
% (3) Type - 'Gaussian', or 'Butterworth'
% (4) Transform Type 'DCT', 'FFT'
% (5) If Butterworth is selected, must select order

% Outputs:
% (1) Filtered Image 
% (2) Centered Kernel 
% (3) Filtered Image (Freq domain)
% (4) Original Image in Freq domain
%
%Written by Nikeet Pandit

function [Im_Filter, Im_Diff, Im_Filter_Freq, Im_DFT_return, Pt, Telapsed] = lowpass_im(Im, D0, kernel_type, transform_type, varargin)

Tstart = tic;

transform_type = lower(transform_type); 
kernel_type = lower(kernel_type); 

if length(varargin) == 1 %Setting order for Butterworth if selected
    n = varargin{1};
else
    n = []; 
end

%% Determine image size and calculate padded image size 
[M, N] = size(Im); P = M*2; Q = N*2; 

%% Construct kernel (frequency domain) to be used for filtering
H = kernel_construct(D0 , P, Q, kernel_type, n); %function is appended below

%% Take DFT of the image with padding specified by P and Q

if transform_type == "fft"
    Im_DFT = fft2((Im),P,Q); 
    Im_DFT_return = abs(fftshift((Im_DFT))); %centered spectrum 
elseif transform_type == "dct"
    Im_DFT = dct2(Im,P,Q); 
    Im_DFT_return = fftshift(Im_DFT); 
else
    print("Improper Selection");
end
   

%% Filter in frequency domain 
Im_Filter = (H.*Im_DFT);

if transform_type == "fft"
    Im_Filter_Freq = fftshift(abs(Im_Filter)); %returning filtered spectrum 
elseif transform_type == "dct"
    Im_Filter_Freq = fftshift((Im_Filter)); %returning filtered spectrum 
end

%% Check Power
Pt = (sum(abs(Im_Filter),'all')/(sum(abs(Im_DFT),'all')))*100; 

%%Isolate real components only
if transform_type == "fft"
    Im_Filter = (real(ifft2(Im_Filter))); 
elseif transform_type == "dct"
    Im_Filter = uint8(idct2(Im_Filter));
end
   
%% Crop Image to remove padding 
Im_Filter = (Im_Filter(1:M, 1:N)); %returning cropped filtered image
H = fftshift(H);  %returning centered (and cropped) kernel function 

Telapsed = toc(Tstart); 

Im_Diff = Im - Im_Filter; %High pass filter of original and filtered

end
