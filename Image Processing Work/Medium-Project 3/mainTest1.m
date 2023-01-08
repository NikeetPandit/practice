%Assignment 3 by Nikeet Pandit
%GS/MATH 6920 Harmonic Analysis and Image Processing
%

%% Test 1) 
%%
% In test 1), we will look at periodic noise reduction using frequency
% domain filtering. More specifically, we will look at notch filtering in
% the frequency domain. First, the image will be read in and periodic noise
% is added to the image. Periodic noise is added to the image in the
% frequency domain for simplicity. Thereafter, a notch reject filter is
% constructed to reject the added periodic noise. Codes for constructing
% the notch reject kernels have been written for this test 1). Metrics to
% test the efficacy of the notch filter - mean-square error and structural similarity index are
% employed to compare the original and filtered image. Other types of noise are tested to see the resultant spectrum.
% 
% When the noise is not periodic, it is clear that the image could likely
% be improved with a lowpass filter to smooth out high frequency
% variations. However, other spatial algorithms (for example adaptive
% spatial filtering likely may perform better) and preserve more edge
% details that the low-pass frequency domain filtering is not capable of.
% More so, the effect of different sorts of noise in the spectrum is
% evident in the visual change of the noise added image spectrum compared
% against the original, however, it is impossible to exactly detect where
% in the spectrum the noise is exactly contained because it is not exactly
% periodic at one frequency where it can easily be removed. 
% When the noise is periodic, it can easily be removed with some
% investigation of the spectrum as will be seen. 
%  
% 

clearvars; close all
%% Reading Image 
Im = imread("rose512.tif"); 

%% Adding Periodic Noise

%--- FFT of image and centering
Im_DFT = fftshift(fft2((double(Im)))); 

%--- Calculating distance matrix from centre of frequency domain rectangle
[M,N] = size(Im); 
u = 0:M-1; %frequency components
v = 0:N-1; 
[U, V] = meshgrid(u,v); 
D = hypot(U-M/2, V-N/2); %calculating distances

%--- Extracting top left and bottom right index from distance matrix equal to 348
[row,col] = find(D == 348); 
TopLeft = [row(1), col(1)]; %TopLeft and BotRight and on opposite corners of D matrix
BotRight = [row(end), col(end)]; 

%--- Adding noise at these locations 
Noise = 2e7; 
Im_DFT(TopLeft(1),TopLeft(2)) = Im_DFT(TopLeft(1),TopLeft(2)) + Noise; 
Im_DFT(BotRight(1),BotRight(2)) = Im_DFT(BotRight(1),BotRight(2)) + Noise; 

%--- Inverting by IFFT2
Im_Noisy = ifft2(fftshift(Im_DFT)); 

%--- Plotting original and noise added spectrum
fh = figure(1); fh.WindowState = 'maximized';
subplot_tight(1,2,1)
imagesc((abs(real(fftshift(fft2((double(Im))))))));
colorbar; title('Original Spectrum'); 
subplot_tight(1,2,2); imagesc((abs(real(Im_DFT)))); colorbar; title('Noise Added Spectrum'); 

%%
%--- Setting cut off frequency and location of notch filter offset 
%offset is with respect to centre of frequency centre
D0 = 5; 
Vk = 3; Uk = 16; %determined as function of where periodic noise and visually from spectrum
kernel = kernel_construct_notch(D0, Uk, Vk, M, N, 'gaussian'); 
fh = figure(2); fh.WindowState = 'maximized';
h = surf(fftshift(kernel)); title('Notch Filter Gaussian Kernel');

%--- Filtering Periodic Noisy Image with Notch Filter and recovering image %in spacetial domain
Im_Filtered = ifft2(fftshift((fftshift(fft2(Im_Noisy)).*kernel)));

%--- Plotting All Results
fh = figure(3); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(uint8(Im)); title('Original Image');
subplot_tight(1,3,2); imshow(uint8(Im_Noisy)); title('Periodic Noise Added'); 
subplot_tight(1,3,3); imshow(uint8(Im_Filtered)); title('Notch Filtered'); 

%%
%--- Testing Notch Filter Periodic Noise Removal in Frequency Domain with Other noise sources
fh = figure(4); fh.WindowState = 'maximized';
Im_gaus = imnoise(Im,'gaussian',0,0.1); %mean 0 var 0.1
Im_sp = imnoise(Im,'salt & pepper',0.05); %salt and pepper noise density 0.05
Im_uform = double(Im)+10*rand(M,N); %uniform noise a = 0, b = 10

subplot_tight(2,2,1); imshow(uint8(Im)); title('Original Image');
subplot_tight(2,2,2); imshow(uint8(Im_gaus)); title('Gaussian Noise Added'); 
subplot_tight(2,2,3); imshow(uint8(Im_sp)); title('Salt & Pepper'); 
subplot_tight(2,2,4); imshow(uint8(Im_uform)); title('Uniform Noise'); 

%%
%--- Plotting Amp Spectrums
fh = figure(5); fh.WindowState = 'maximized';
%--- Calc. amplitude spec in dB
ImSpec = 20*log10(abs(real(fftshift(fft2(Im))))); %isolating amplitude component and taking abs value
Im_gausSpec = 20*log10(abs(real(fftshift(fft2(Im_gaus))))); 
Im_spSpec = 20*log10(abs(real(fftshift(fft2(Im_sp)))));
Im_uformSpec = 20*log10(abs(real(fftshift(fft2(Im_uform))))); 

subplot_tight(2,2,1); imagesc(ImSpec); colorbar; title('Original Image Spec');
subplot_tight(2,2,2); imagesc(Im_gausSpec); colorbar; title('Gaussian Noise Added Spec'); 
subplot_tight(2,2,3); imagesc(Im_spSpec); colorbar; title('Salt & Pepper Spec'); 
subplot_tight(2,2,4); imagesc(Im_uformSpec);  colorbar; title('Uniform Noise Spec'); 

%%
% 
% As we have seen, the periodic noise was added to the image which in the
% spatial domain the effect of this added noise completely corrupted the
% image. Comparing the spectrum with the original image and the periodic
% noise added image, the added periodic noise manifests itself in
% spectrally as two bursts of high amplitude energy. Thereafter, a notch
% filter can be constructed whose centre frequency is merely shifted to the
% location of these high frequency bursts. Some inspection of the spectrum
% was required to see where to shift the frequencies of the notch filter
% because the shift Uk and Vk is horizontal and vertical shifting with
% respect to the centre of the frequency rectangle. A gaussian kernel was
% used with a very sharp cut off because it is smooth, so no ringing
% artefacts are introduced in the filtering process.
% 
% The MSE for the corrupted image (without filtering) is 11642 which is
% very high. The filtered MSE has a result of 4.5 which is significantly
% less by a factor of 2500. 
% The structuraly similar index, another measure of image quality, shows a
% value of 0.8 (closer to 1 is better) for the filtered result and a value
% of 0 for the corrupted image. Clearly, the effect of the removal of
% periodic noise is very efficient using notch filters. 

% When other types of noise are added, namely salt and pepper, gaussian,
% and uniform, there is a change in the resultant spectrum compared to the
% original, as showed by the introduction of additional higher frequency
% component amplitude power. Although the uniform noise spectrum looks very
% similar to the original. This example to show is that that while a low
% pass filter would likely improve image quality and may remove some of the
% added noise (at the cost of edge detail), the effect of the noise that is
% not completely periodic is not as obvious to spot in the spectrum and
% therefore not easy (or impossible) to be completely removed using
% frequency domain filtering alone. In these cases, perhaps spatial filters
% such as the adaptive mean filters would be a better option. Choosing
% variable amounts of noise yields similar results as the test shown in
% this test.



