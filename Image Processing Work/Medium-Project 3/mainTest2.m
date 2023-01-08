%Assignment 3 by Nikeet Pandit
%GS/MATH 6920 Harmonic Analysis and Image Processing
%
clearvars; close all
%%
% 
% In this test, two degredation (or point spread functions) will be tested.
% Namely, blurring function and a motion function. For both the blurring
% function and the motion function, there will be tests with salt and
% pepper noise, gaussian noise, and uniform noise. For the motion
% degradation function, there will be more noise inputted then the blurring
% function. Differences will be examined. Estimates for the noise function
% will be made and compared to the actual noise inputted. The filters
% tested for this series of experiements will be the inverse filter,
% constrained LS filter, and the weiner filter. 
% 


%% Reading Image 
Im = im2double(imread("rose512.tif")); 

%% Simulate to degredation functions 
%--- Blurring degredation function 
hsize = 10; sigma = 5; 
H_blur = fspecial('gaussian', hsize, sigma);

%--- Applying blur degredation function to image
Im_blur = imfilter(Im, H_blur,'conv','circular'); 

%--- Motion degredation function
len = 25; theta =  60; 
H_motion = fspecial('motion', len, theta);

%--- Applying motion degredation function to image
Im_motion = imfilter(Im, H_motion,'conv','circular');

%--- Plotting degraded images alongside original 
fh = figure(1); fh.WindowState = 'maximized';; subplot_tight(1,3,1); imshow(Im); title('Original'); 
subplot_tight(1,3,2); imshow(Im_blur); title('Blurred'); 
subplot_tight(1,3,3); imshow(Im_motion); title('Motion'); 
%%
%% Add noise to degraded images 
% S&P, Gaussian, Uniform Noise are added
ab_gauss = [0, 0.01]; %mean 0, variance 0.01 gaus noise
d_sp = 0.05; %density 0.05 S&P noise (equal probability); 
ab_uform = [0, 0.1]; %a = 0, b = 0.1 uniform noise

%--- adding noise to blurred degraded image
[~, Im_SP_blur, Im_Gaus_blur, Im_Uform_blur] = add_noise_func(Im_blur, ab_gauss, d_sp, ab_uform);

fh = figure(2); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_SP_blur); title('S&P Noise + Blur'); 
subplot_tight(1,3,2); imshow(Im_Gaus_blur); title('Gaussian Noise + Blur'); 
subplot_tight(1,3,3); imshow(Im_Uform_blur); title('Uniform Noise + Blur'); 
%%

%--- S&P, Gaussian, Uniform Noise are added
%--- Please note there is much more noise added for the motion PSF then blur PSF
ab_gauss1 = [0, 0.01]*10; %mean 0, variance 0.1 gaus noise
d_sp1 = 0.1; %density 0.1 S&P noise (equal probability); 
ab_uform1 = [0, 0.15]; %a = 0, b = 0.15 uniform noise

%--- adding noise to motion degraded image
[~, Im_SP_motion, Im_Gaus_motion, Im_Uform_motion] = add_noise_func(Im_motion, ab_gauss1, d_sp1, ab_uform1);

fh = figure(3); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_SP_motion); title('S&P Noise + Motion'); 
subplot_tight(1,3,2); imshow(Im_Gaus_motion); title('Gaussian Noise + Motion'); 
subplot_tight(1,3,3); imshow(Im_Uform_motion); title('Uniform Noise + Motion'); 
%%
%% Estimate noise variance from degraded (and noise added) image

%--- I use a low pass filter on the degraded / noisy images
% After, I use a difference filter with the original degraded (with noise)
% imaged and the low pass filter to highlight high frequency differences
% between the two. 

% Since degredation functions are introducing blurs we would expect any
% high frequency difference between the degredated image and the smoothed
% degredation image would attributed to the noise. Therefore, I use this
% as noise estimate and compared to the actual value of the noise variances to
% see its efficacy as a noise variance estimate. 

%--- Smoothing blur degredation image with noise added
%--- And returning difference filter (high pass)
%--- ... = original(degraded  %noise) - smoothed(degraded + noise)

DO = 120; %cut off frequency
[Im_Filter_Gaus_blur, Im_Diff_Gaus_blur, ~, ~, ~, ~] = lowpass_im(Im_Gaus_blur, DO, 'gaussian', 'FFT'); %Smoothing gaussian noise added image
[Im_Filter_SP_blur, Im_Diff_SP_blur, ~, ~, ~, ~] = lowpass_im(Im_SP_blur, DO, 'gaussian', 'FFT'); %Smoothing S&P noise added image
[Im_Filter_Uform_blur, Im_Diff_Uform_blur, ~, ~, ~, ~] = lowpass_im(Im_Uform_blur, DO, 'gaussian', 'FFT'); %Smoothing uniform noise added image

%--- Extracting quadrants of image (lower right, lower left, upper right, upper left)
%... then averaging to take noise estimate where there is only background

Im_Diff_Gaus_blur = quadrant_avg_func(Im_Diff_Gaus_blur);
Im_Diff_Uform_blur = quadrant_avg_func(Im_Diff_Uform_blur); 

%--- For S&P taking only upper left quadrant 
Im_Diff_SP_blur = rescale(uint8((Im_Diff_SP_blur(1:100,1:100))),0,255);

%--- Plotting Noise Patches Where Noise Variance is Est from
fh = figure(4); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_Diff_SP_blur); title('Salt Noise + Blur');
subplot_tight(1,3,2); imshow(Im_Diff_Gaus_blur); title('Gaussian Noise + Blur'); 
subplot_tight(1,3,3); imshow(Im_Diff_Uform_blur); title('Uniform Noise + Blur'); 
suptitle('Noise Patch Plots (noise var is est. from here)')
%%
%--- Calculating Pepper Noise (probability) all other intensities except 255 set to 0
ProbOfSalt = numel(find(Im_Diff_SP_blur == 255))/numel(Im_Diff_SP_blur); 

%--- Equal salt and pepper noise assumption 
d_sp_est = ProbOfSalt *2; 

%--- Estimating noise variance from difference filtered image
var_est_gaus_blur = var(Im_Diff_Gaus_blur(:)); 
var_est_SP_blur = SP_variance_calc(d_sp_est);
var_est_Uform_blur = var(Im_Diff_Uform_blur(:)); 

Est_Blur_Var = [var_est_gaus_blur var_est_SP_blur var_est_Uform_blur]; 

%--- Getting real noise variance based on inputted values
var_real_gaus_blur = ab_gauss(2);
var_real_SP_blur =   SP_variance_calc(d_sp); 
var_real_Uform_blur = uform_var_calc(ab_uform(1),ab_uform(2)); 

Real_Blur_Var = [var_real_gaus_blur var_real_SP_blur var_real_Uform_blur];

%--- Calculating the difference between estimate and actual
Diff_Blur_Var = abs(Real_Blur_Var - Est_Blur_Var); 

%% Estimate noise variance from degraded (and motion added) image

%--- Smoothing blur degredation image with noise added
%--- And returning difference filter (high pass)
%--- ... = original(degraded  %noise) - smoothed(degraded + noise)

DO = 120; %cut off frequency
[Im_Filter_Gaus_motion, Im_Diff_Gaus_motion, ~, ~, ~, ~] = lowpass_im(Im_Gaus_motion, DO, 'gaussian', 'FFT'); %Smoothing gaussian noise added image
[Im_Filter_SP_motion, Im_Diff_SP_motion, ~, ~, ~, ~] = lowpass_im(Im_SP_motion, DO, 'gaussian', 'FFT'); %Smoothing S&P noise added image
[Im_Filter_Uform_motion, Im_Diff_Uform_motion, ~, ~, ~, ~] = lowpass_im(Im_Uform_motion, DO, 'gaussian', 'FFT'); %Smoothing uniform noise added image

%--- Extracting quadrants of image (lower right, lower left, upper right, upper left)
%... then averaging to take noise estimate where there is only background

Im_Diff_Gaus_motion = quadrant_avg_func(Im_Diff_Gaus_motion);
Im_Diff_Uform_motion = quadrant_avg_func(Im_Diff_Uform_motion); 

%--- For S&P taking only upper left quadrant 
Im_Diff_SP_motion = rescale(uint8((Im_Diff_SP_motion(1:100,1:100))),0,255);

%--- Plotting Noise Patches Where Noise Variance is Est from
fh = figure(5); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_Diff_SP_motion); title('S&P Noise + Motion');
subplot_tight(1,3,2); imshow(Im_Diff_Gaus_motion); title('Gaussian Noise + Motion'); 
subplot_tight(1,3,3); imshow(Im_Diff_Uform_motion); title('Uniform Noise + Motion'); 
suptitle('Noise Patch Plots (noise var is est. from here)')
%%
%--- Calculating Pepper Noise (probability) all other intensities except 255 set to 0
ProbOfSalt = numel(find(Im_Diff_SP_motion == 255))/numel(Im_Diff_SP_motion); 

%--- Equal salt and pepper noise assumption 
d_sp_est = ProbOfSalt *2; 

%--- Estimating noise variance from difference filtered image
var_est_gaus_motion = var(Im_Diff_Gaus_motion(:)); 
var_est_SP_motion = SP_variance_calc(d_sp_est);
var_est_Uform_motion = var(Im_Diff_Uform_motion(:)); 

Est_motion_Var = [var_est_gaus_motion var_est_SP_motion var_est_Uform_motion]; 

%--- Getting real noise variance based on inputted values
var_real_gaus_motion = ab_gauss1(2);
var_real_SP_motion =   SP_variance_calc(d_sp1); 
var_real_Uform_motion = uform_var_calc(ab_uform1(1),ab_uform1(2)); 

Real_motion_Var = [var_real_gaus_motion var_real_SP_motion var_real_Uform_motion]; 

%--- Calculating the difference between estimate and actual
Diff_motion_Var = abs(Est_motion_Var-Real_motion_Var); 

%% Testing Inverse Filter 
% Inverse filter uses the image and the PSF
% No estimated noise is provided, except knowledge of the PSF is required. 
% Tests will be done for the two PSF with noise and without noise to see
% when the weiner filter works well. 

%--- Test without noise and use wiener filter 
Im_recov_Blur = deconvwnr(Im_blur, H_blur); 
Im_recov_Motion = deconvwnr(Im_motion, H_motion);

%--- Test without noise but adding noise (uncertainty) to knowledge of PSF 
H_blur_noisy = H_blur + 0.0005*rand(size(H_blur,1),size(H_blur,2));
H_motion_noisy = H_motion .*(1e-1*rand(size(H_motion,1),size(H_motion,2))) + H_motion;

Im_recov_Blur_uncertain = deconvwnr(Im_blur, H_blur_noisy);
Im_recov_Motion_uncertain = deconvwnr(Im_motion, H_motion_noisy);

fh = figure(6); fh.WindowState = 'maximized';
subplot_tight(2,2,1); imshow(Im_recov_Blur); title('Blurred PSF (known exact)'); 
subplot_tight(2,2,2); imshow(Im_recov_Motion); title('Motion PSF (known exact)'); 
subplot_tight(2,2,3); imshow(Im_recov_Blur_uncertain); title('Blurred PSF (uncertain)'); 
subplot_tight(2,2,4); imshow(Im_recov_Motion_uncertain); title('PSF (uncertain)'); 
suptitle('Restorations'); 
%%

% %
%  As we can see, exact knowledge of the PSF is essential. Even a bit of
%  noise added to the PSF causes distortion to restoration of the original
%  image for inverse filtering. When knowledge of the PSF is known exactly,
%  inverse filtering provides a MSE of 6.8e-7 and a SSI (structural similarity
%  index) of 0.99, showing a near perfect score. When visually comparing these
%  restored with the original image there is no visual difference either.
%  So, perfect reconstruction is achieved with perfect knowledge of the
%  PSF under the prescence of no noise. When, however, the PSF are not
%  known exactly the MSE increases to ~5.5e-4 for both PSF
%  and the SSI decreases to 0.93 for the blur PSF and 0.84 for the motion
%  PSF. Visually, there is also lots of distortion present. 
% 

%--- Test with noise for blur PSF
Im_recov_Blur_SP = deconvwnr(Im_SP_blur, H_blur); 
Im_recov_Blur_Gaus = deconvwnr(Im_Gaus_blur, H_blur); 
Im_recov_Blur_Uform = deconvwnr(Im_Uform_blur, H_blur); 

fh = figure(7); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_recov_Blur_SP); title('S&P Noise added'); 
subplot_tight(1,3,2); imshow(Im_recov_Blur_Gaus); title('Gaus Noise added'); 
subplot_tight(1,3,3); imshow(Im_recov_Blur_Uform); title('Uniform Noise added'); 
suptitle('Restorations From Blur PSF'); 
%%

%% 
%  As we can see, inverse filtering completely breaks down under the
%  prescence of noise and the image cannot be restored at all! The same
%  process happens when using the motion PSF and those results will not be
%  appended for brevity. Overall, these experiments demonstrate the linear
%  filtering is great when the PSF is known exactly under the prescence of
%  no noise to restore the original image. 

%% Testing Wiener Filtering
 
%--- Test with noise for blur PSF

Im_recov_Blur_SP = deconvwnr(Im_SP_blur, H_blur,var_real_SP_blur/var(Im(:))); 
Im_recov_Blur_Gaus = deconvwnr(Im_Gaus_blur, H_blur, var_real_gaus_blur/var(Im(:))); 
Im_recov_Blur_Uform = deconvwnr(Im_Uform_blur, H_blur, var_real_Uform_blur/var(Im(:))); 

fh = figure(8); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_recov_Blur_SP); title('S&P Noise added'); 
subplot_tight(1,3,2); imshow(Im_recov_Blur_Gaus); title('Gaus Noise added'); 
subplot_tight(1,3,3); imshow(Im_recov_Blur_Uform); title('Uniform Noise added'); 
suptitle('Wiener - Blur PSF with known NSR'); 
%%

%--- Test with noise for blur PSF using estimated SNR
Im_recov_Blur_SP = deconvwnr(Im_SP_blur, H_blur,var_est_SP_blur/var(Im_SP_blur(:))); 
Im_recov_Blur_Gaus = deconvwnr(Im_Gaus_blur, H_blur, var_est_gaus_blur/var(Im_Gaus_blur(:))); 
Im_recov_Blur_Uform = deconvwnr(Im_Uform_blur, H_blur, var_est_Uform_blur/var(Im_Uform_blur(:))); 

fh = figure(9); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_recov_Blur_SP); title('S&P Noise added'); 
subplot_tight(1,3,2); imshow(Im_recov_Blur_Gaus); title('Gaus Noise added'); 
subplot_tight(1,3,3); imshow(Im_recov_Blur_Uform); title('Uniform Noise added'); 
suptitle('Wiener - Blur PSF with est.  noise power and signal power'); 
%%
%%
% As these tests show, for Wiener filtering, exact knowledge of the NSR,
% that is to say, knowledge of the variance of the original signal and its
% noise is critical. Although the estimated noise is very close to the
% actual inputted noise (estimate noise variance has less then 0.01 
% absolute difference from the real or even less then for all estimates). 
% Even with knowledge of the input signal variance (original image) and
% variances for all the noise (i.e., inputted parameters), the
% reconstruction is not perfect. For S&P noise the SSI is 0.72 and the MSE
% is 0.02. For Gaussian noise, the SSI is 0.65 and the MSE is 0.0035. For
% unform noise, the SSI is 0.3 and the MSE is 0.0054. Interestingly, there
% is an inverse relationship for the tests between the SSI and MSE.
% Theoretically, the for a better image the SSI should be higher the MSE
% should be lower, however, this relationship is not maintained as
% expected. Visually, however, the images are not very well reconstructed
% and the effect of noise is very obvious. For the test using estimated NSR
% ratio, the effect is visual as well as reflected strongly in the SSI with
% much lower values than compared to when using the known NSR ratio. 
% 


% --- Try improving the weiner estimate by using acf of noise and original image
% ifft2 of power spectrum yields acvf... normalized by rms squared then gives acf

%--- Noise is determined by subtracting noise added image by original
acf_Image = (ifft2((fft2(abs(Im).^2)/numel(Im))))/(rms(Im(:))^2);
acf_NoiseSP = (ifft2(fft2(abs(Im_SP_blur - Im).^2))/numel(Im))/rms(Im(:))^2;
acf_NoiseGaus = (ifft2(fft2(abs(Im_Gaus_blur - Im).^2))/numel(Im))/rms(Im(:))^2;
acf_NoiseUform = (ifft2(fft2(abs(Im_Uform_blur - Im).^2))/numel(Im))/rms(Im(:))^2;

%--- Using Wiener filter
Im_recov_Blur_SP = deconvwnr(Im_SP_blur, H_blur,acf_NoiseSP, acf_Image); 
Im_recov_Blur_Gaus = deconvwnr(Im_Gaus_blur, H_blur,acf_NoiseGaus, acf_Image); 
Im_recov_Blur_Uform = deconvwnr(Im_Uform_blur, H_blur,acf_NoiseUform, acf_Image); 

fh = figure(10); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_recov_Blur_SP); title('S&P Noise added'); 
subplot_tight(1,3,2); imshow(Im_recov_Blur_Gaus); title('Gaus Noise added'); 
subplot_tight(1,3,3); imshow(Im_recov_Blur_Uform); title('Uniform Noise added'); 
suptitle('Wiener Restoration - Blur PSF with ACF of noise and signal'); 
%%

%%
% When inputting the ACF (estimated through the IFFT2) visually the results
% faired much better, especially for the uniform noise. 
% Now, using only the SSI for the performance metric
% (seemed to fair better and correspond better my visual preference than
% MSE), the SSI for the S&P noise, gaussian noise, and uniform noise are
% 0.68, 0.61, 0.6. Although not reflected in the SSI, I prefer the the
% restoration using the ACF for gaussian noise and the S&P noise as well.
% The uniform noise restoration does particularly well compared to all
% others. 

%% Testing Wiener Filtering with Motion PSF (and more noise)
% --- Tests will only be done with the ACF for the signal and noise which
% seemed to fair the best

%--- Noise is determined by subtracting noise added image by original
acf_Image = (ifft2((fft2(abs(Im).^2)/numel(Im))))/(rms(Im(:))^2);
acf_NoiseSP = (ifft2(fft2(abs(Im_SP_motion - Im).^2))/numel(Im))/rms(Im(:))^2;
acf_NoiseGaus = (ifft2(fft2(abs(Im_Gaus_motion - Im).^2))/numel(Im))/rms(Im(:))^2;
acf_NoiseUform = (ifft2(fft2(abs(Im_Uform_motion - Im).^2))/numel(Im))/rms(Im(:))^2;

%--- Using Wiener filter
Im_recov_motion_SP = deconvwnr(Im_SP_motion, H_motion,acf_NoiseSP, acf_Image); 
Im_recov_motion_Gaus = deconvwnr(Im_Gaus_motion, H_motion,acf_NoiseGaus, acf_Image); 
Im_recov_motion_Uform = deconvwnr(Im_Uform_motion, H_motion,acf_NoiseUform, acf_Image); 

fh = figure(11); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_recov_motion_SP); title('S&P Noise added'); 
subplot_tight(1,3,2); imshow(Im_recov_motion_Gaus); title('Gaus Noise added'); 
subplot_tight(1,3,3); imshow(Im_recov_motion_Uform); title('Uniform Noise added'); 
suptitle('Wiener Restoration - Motion PSF with ACF of noise and signal'); 
%%
%%
% In this case, although the noise removal visually quite good and
% and comparable to the previous example (which has less noise)
% the effect the motion is not completely removed and is present even with
% the full knowledge of the motion PSF, despite it being able to be removed
% in inverse filtering (with no noise present). All of this means and goes
% to show that for Wiener filtering knowledge of the signal and noise ACF
% and knowledge of the PSF is essential to optimally restore the image. 

%% Testing Least Squares Constrained Filtering 
% Since for Wiener filtering the ACF for the noise and signal must be
% known, this brings additional difficulty as these may not be known and
% estimating them can be quite difficult. As seen in the above examples,
% estimates that differ even slightly from the "true" values yield restoration results
% that are not very effective. In the CLS, tests will be focused to see the
% effect of noise and see how it compares against the Wiener filtering
% restoration result and how sensitive the CLS is to noise power estimate
% uncertainty, as well as no knowledge of noise power. 

%--- Using LSC filter with no input for noise power (Blur PSF)
Im_recov_motion_SP = deconvreg(Im_SP_blur, H_blur); 
Im_recov_motion_Gaus = deconvreg(Im_Gaus_blur, H_blur); 
Im_recov_motion_Uform = deconvreg(Im_Uform_blur, H_blur); 

fh = figure(12); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_recov_motion_SP); title('S&P Noise added'); 
subplot_tight(1,3,2); imshow(Im_recov_motion_Gaus); title('Gaus Noise added'); 
subplot_tight(1,3,3); imshow(Im_recov_motion_Uform); title('Uniform Noise added'); 
suptitle('CLS Restoration - Blur PSF w/o input for noise power'); 
%%
%%
%Comparing inverse filtering with CLS, where both filters require only the
%input image and knowledge of the PSF, the results for the CLS are
%significantly better as the rose structure can at least be uncovered all
%the noise, compared to absolutely no rose structure in the inverse
%filtering example. 

%--- Using LSC filter with no input for noise power (Motion PSF) 
Im_recov_motion_SP = deconvreg(Im_SP_motion, H_motion); 
Im_recov_motion_Gaus = deconvreg(Im_Gaus_motion, H_motion); 
Im_recov_motion_Uform = deconvreg(Im_Uform_motion, H_motion); 


fh = figure(13); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_recov_motion_SP); title('S&P Noise added'); 
subplot_tight(1,3,2); imshow(Im_recov_motion_Gaus); title('Gaus Noise added'); 
subplot_tight(1,3,3); imshow(Im_recov_motion_Uform); title('Uniform Noise added'); 
suptitle('CLS Restoration - Motion PSF w/o input for noise power'); 
%%
%%
%With more noise added and the PSF function set to motion the rose
%structure cannot be uncovered without any noise power information. In the
%uniform noise example, the structure can just barely be distinguished.
%Nevertheless, the CLS is a significant improvement over the inverse
%example, despite requiring the same number of input parameters. 

%--- Using LSC filter with using "true" noise power values (Blur PSF)
n = numel(Im); 
Im_recov_motion_SP = deconvreg(Im_SP_blur, H_blur,var_real_SP_blur*n); 
Im_recov_motion_Gaus = deconvreg(Im_Gaus_blur, H_blur, var_real_gaus_blur*n); 
Im_recov_motion_Uform = deconvreg(Im_Uform_blur, H_blur, var_real_Uform_blur*n); 

fh = figure(14); fh.WindowState = 'maximized';
subplot_tight(1,3,1); imshow(Im_recov_motion_SP); title('S&P Noise added'); 
subplot_tight(1,3,2); imshow(Im_recov_motion_Gaus); title('Gaus Noise added'); 
subplot_tight(1,3,3); imshow(Im_recov_motion_Uform); title('Uniform Noise added'); 
suptitle('CLS Restoration - Blur PSF w/ input for noise power (real noise)'); 
%%



%%
% This example is showing clearly there is some issue with the variance
% calculation for the variance of S&P noise. Even though I double checked
% the formula, I am not sure how to debug it. Including noise power into
% the image brings a significant improvement, however, the blurring effect
% is quite evident. Comparing visually, the effect is similar for uniform
% noise to the Wiener filter (except the Weiner filter requires ACF of
% signal as well, while CLS only requires np). While the Gaussian noise
% image seems to have the noise completely removed, all the high frequency
% information is obviously very smoothed away. Now, there will be some
% experimentation with the value of noise power to yield a sharpened
% result. Since there is some error with the S&P variance calculation, only
% Gaussian and Uniform noise will be tested now. 

fraction = 0.65; %only include fraction of noise power to sharpen result
fraction1 = 0.995;
Im_recov_motion_Gaus = deconvreg(Im_Gaus_blur, H_blur, (var_real_gaus_blur*n)*fraction); 
Im_recov_motion_Uform = deconvreg(Im_Uform_blur, H_blur, fraction1*(var_real_Uform_blur*n)); 

figure(15); 
subplot_tight(1,2,1); imshow(Im_recov_motion_Gaus); title('Gaus Noise added'); 
subplot_tight(1,2,2); imshow(Im_recov_motion_Uform); title('Uniform Noise added'); 
suptitle('CLS Restoration - Blur PSF w/ input for noise power (adjusted noise)'); 
%%

%%
% As we can see, variably adjusting the noise power can yield a better
% result. This is especially true for the Gaussian noise added example. For
% the uniform noise example, there is little to no improvement with
% variably changing the noise power. This shows that absolute knowledge of
% the noise power is not as necessary as the Wiener filter, and is thus
% more flexible and can be parameterized to yield better results. A
% also, by being able to adjust the noise power level to
% trade blurring for sharpening and vice versa, the filter offers great
% flexibility. 

%--- LSC Filter with noise power values (Motion PSF) 
%--- Using LSC filter with using "true" noise power values (Blur PSF)
n = numel(Im); 
Im_recov_motion_Gaus = deconvreg(Im_Gaus_motion, H_motion, (var_real_gaus_motion*n)*0.5); 
Im_recov_motion_Uform = deconvreg(Im_Uform_motion, H_motion, var_real_Uform_motion*n); 

figure(16); 
subplot_tight(1,2,1); imshow(Im_recov_motion_Gaus); title('Gaus Noise added'); 
subplot_tight(1,2,2); imshow(Im_recov_motion_Uform); title('Uniform Noise added'); 
suptitle('CLS Restoration - Motion PSF w/ input for noise power'); 
%%
%%
% For both the Wiener filter and the LSC filter, when the degredation
% filter is motion, that is to see there is a blurring effect on top of a
% the movement distrotion, the restoration effect of the image is not as
% effective in both the removal of the noise and removal of the
% degredadation function (deconvolution). This is to say, that the more
% complicating the degredation function, the more difficult it is for these
% filtering methods to restore the original image. Visually, the effect of
% restoration effect of the CLS for the uniform noise introduced additional
% distortion and ripppling throughout the image, which is very apparent in
% the background. The Gaussian noise added example also was very blurred
% and even with tinkering the noise power value, would result in either an
% image that is too nosiy or too blured. All of this to say is that the
% Wiener filter is a much better option provided a good estimate for signal
% ACF and the noise ACF can be provided. This can be verified visaully by
% comparing figures (11) and figures (16) and also improvements in the SSI
% for the Wiener constructed images. However, when these results are not
% known, or the NSR is not very well estimated, the CLS offers a much
% better construction result. 










