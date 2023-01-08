%Assignment 2 by Nikeet Pandit
%GS/MATH 6920 Harmonic Analysis and Image Processing

% --- Question List --- %
% 4.7a) Read in "blurry_moon" and sharpen using unsharp masking in 
%       FREQUENCY DOMAIN
%
% 4.7b) Improve sharpness using highboost filtering
%
% --- Question List --- %

%% 4.7 (a) and (b)
close all; clearvars; addpath('functions');

%% Read the image
Im = imread("blurry_moon.tif"); 

%% Determine image size and calculate padded image size 
[M, N] = size(Im); P = M*2; Q = N*2; 

%% Construct kernel (frequency domain) to be used for filtering
D0 = 10; %cut off frequency
H = kernel_construct(D0,P,Q,'Gaussian'); %function is uploaded to GitHub

%% Take DFT of the image with padding specified by P and Q
Im_DFT = fft2((Im),P,Q); 

%% Filter in frequency domain and isolate real components only
Im_Filter = real(ifft2(H.*Im_DFT)); 

%% Crop Image to remove padding 
Im_Filter = Im_Filter(1:M, 1:N); 

%% Frequency Domain Unsharp Masking (a)
mask = double(Im) - Im_Filter; 
k = 1; 
G_unsharp = (1 + k*(1-H)).*Im_DFT;  %Expression in frequency domain
g_unsharp = real(ifft2(G_unsharp)); %Converting back to spatial domain

%% Frequency Domain Highboost filtering
k = 1.5; 
G_highboost = (1 + k*(1-H)).*Im_DFT; %Expression in frequency domain
g_highboost = real(ifft2(G_highboost)); %Converting back to spatial domain

%% Plotting
figure(1); 
subplot(2,2,1); 
h = surf(fftshift(H)); title('Centred Gaussian Transfer Func'); 
set(h,'LineStyle','none'); 
subplot(2,2,2); 
imshow(Im); title('Original Image'); 
subplot(2,2,3); 
imshow(uint8(Im_Filter)); title('Blurred Image'); 
subplot(2,2,4); 
imshow(uint8(mask)); title('Filter Mask'); 

figure(2); 
subplot(1,2,1); imshow(uint8(g_unsharp(1:M,1:N)),[]); title('Unsharp Masking'); 
subplot(1,2,2); imshow(uint8(g_highboost(1:M,1:N)),[]); title('Highboost Filter'); 

%% Some Commentary on 4.7
% The effect and improved sharpeness is readily apprent in high-boost
% filtering. Clearly, both unsharp masking and high-boost filtering are
% sharpened with respect to the original image, where the high frequency
% edge components of the images are enhanced. 