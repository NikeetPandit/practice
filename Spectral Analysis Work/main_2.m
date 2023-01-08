%Spectral Analysis Practice and Experimentation 
%
% --- Objective List --- %
% 1) Compare different methods for calculating the PSD.. looking at ACVF,
%FFT, and LSSA.
%
% 2) Investigate how you can have aliasing when the Nyquist condition is
%respected! Known as sub-Nyquist artefacts.
%
% 3) Make the series "gappy" and see the effect it has on the spectrum!
% Observe that LSSA still works! Correpsonds to Q4
%
% 4) More gappy tests. Corresponds to Q5
%
% --- Question List -- %


%% Test 1)
% Compare different methods for power spectrum, PSD
% 1) Will be determined through ACVF / ACF
% 2) Will be determined through square of FFT
% 3) Fourier methods and LSSA methods will be used
%--- No variable freq (already discussed in last week tests)

%---- Constructing series without variable freq, with linear trend
amp = [10 20 420 1000 500]; 
phase = [0 60 20 15 10]; % degrees 
freq_start = [300 430 550 3000 4000]; 
freq_end = freq_start; 
Fs = 10000; 
T = 2;  
SNR_dB = 25;
linear_trend_frac = repmat(1,[1,5]); 
[time, signal] = multiple_synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac);

% ---- Plottign Original Series
figure(1)
plot(time, signal); xlabel('time [s]'); ylabel('volts'); title('Original Synthetic Series'); 

% --- Checking power spectrum (square of FFT, acvf, LSSA) with Linear Trend
figure(2); 
[pxx,w] = periodogram(signal, rectwin(length(signal)), length(signal), Fs,'power');
[freq, amp] = my_FFT(time,signal);
plot(freq,amp)
subplot(1,3,1); plot(w,pxx); title('Spec - FFT'); ylabel('volts^2'); xlabel('Hz'); 
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time, signal,[1, 1, 0, 0]); 
subplot(1,3,2); power = LSSA_helper(CScoeff,'power'); plot(Omega,power); title('Spec - LSSA'); ylabel('volts^2'); xlabel('Hz')
[acf, lags] = autocorr(signal, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(1,3,3); 
[freq, amp1] = my_FFT(time,acvf); 
plot(freq,(amp1*2)); title('Spec - ACVF'); ylabel('volts^2'); xlabel('Hz')

%% Checking power spectrum (square of FFT, acvf, LSSA) without linear trend
close all
% --- Constructing series (same as before) without trend
linear_trend_frac = repmat(0,[1,5]); 
amp = [10 20 420 1000 500]; 
[time, signal] = multiple_synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac);

%--- Plotting Original Series
figure(3)
plot(time, signal); xlabel('time [s]'); ylabel('volts'); title('Original Synthetic Series'); 

%---  Computing Spectrums in dB
figure(4); 
[pxx,w] = periodogram(signal, rectwin(length(signal)), length(signal), Fs,'power');
subplot(1,3,1); plot(w,10*log10(pxx/4)); title('Spec - FFT'); ylabel('dB'); xlabel('Hz')
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time, signal,[1, 1, 0, 0]); 
subplot(1,3,2); power = LSSA_helper(CScoeff,'power'); plot(Omega,10*log10(power/4)); title('Spec - LSSA'); ylabel('dB'); xlabel('Hz')
[acf, lags] = autocorr(signal, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(1,3,3); 
[freq, amp] = my_FFT(time,acvf); 
plot(freq,10*log10((amp*2)/20)); title('Spec - ACVF'); ylabel('dB'); xlabel('Hz')

%---  Computing Spectrums in normal plot
figure(5); 
[pxx,w] = periodogram(signal, rectwin(length(signal)), length(signal), Fs,'power');
subplot(1,3,1); plot(w,pxx); title('Spec - FFT'); ylabel('volts ^2'); xlabel('Hz'); 
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time, signal,[1, 1, 0, 0]); 
subplot(1,3,2); power = LSSA_helper(CScoeff,'power'); plot(Omega,power); title('Spec - LSSA'); ylabel('volts ^2'); xlabel('Hz'); 
[acf, lags] = autocorr(signal, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(1,3,3); 
[freq, amp1] = my_FFT(time,acvf); 
plot(freq,(amp1*2)); title('Spec - ACVF'); ylabel('volts ^2'); xlabel('Hz'); 


% --- Checkign PSD
figure(6); 
[pxx,w] = periodogram(signal, rectwin(length(signal)), length(signal), Fs,'PSD');
subplot(1,3,1); semilogy(w,rdivide(pxx,w)); title('PSD - FFT'); ylabel('volts ^2/Hz'); xlabel('Hz'); 
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time, signal,[1, 1, 0, 0]); 
subplot(1,3,2); power = LSSA_helper(CScoeff,'PSD',Omega); 
semilogy(Omega,power); title('PSD - LSSA'); ylabel('volts ^2/Hz'); xlabel('Hz'); 
[acf, lags] = autocorr(signal, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(1,3,3); 
[freq, amp1] = my_FFT(time,acf); 
semilogy(freq,rdivide((amp1*2),1)); title('PSD - ACF'); ylabel('volts ^2/Hz'); xlabel('Hz'); 

%% Test 2) Sub Nyquist Check for 1 Series

%--- Constructing example case for 1 one series (synethetic)
clearvars; close all; 
amp = 100; 
Fs = 8.25; 
freq_start = 2.5; freq_end = 2.5; 
phase = 0; 
T = 8; 
SNR_dB = 100; 
linear_trend_frac = 0; 
[time, series] = synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac);
Fs = 100000; [time1, series1] = synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac);

%--- Plotting Envolopes based on sub-nyquist rule
Env1 = 1:3:length(series); 
Env2 = 2:3:length(series); 
Env3 = 3:3:length(series); 
figure(1); 
plot(time(Env1), series(Env1),'r'); hold on 
plot(time(Env2), series(Env2),'b'); 
plot(time(Env3), series(Env3),'g'); 
plot(time1,series1,'b');
ylabel('Amplitude'); xlabel('Time [s]'); 
title('Sub Nyquist Effect'); 
legend('Env 1', 'Env 2', 'Env 3', 'g(x)'); 


%% Questions 4
close all; clearvars
% - Series has linear trend, variable freq 
% - Series has also one freq 1/3Fs - 1/4

%--- Constructing Series and Setting Param
amp = [10 20 420 1000 390 500]; 
phase = [0 60 20 15 40 10]; % degrees 
Fs = 10000; 
freq_add = (1/3)*Fs - (1/4); 
freq_start = [300 430 550 3000 freq_add 4500]; 
freq_end = [300 430 750 3000 freq_add 4500]; 
T = 2;  SNR_dB = 100;
linear_trend_frac = repmat(1,[1,6]); 
[time, signal] = multiple_synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac);


%--- Plotting Envolopes based on sub-nyquist rule
figure(1)
Env1 = 1:3:length(signal); 
Env2 = 2:3:length(signal); 
Env3 = 3:3:length(signal); 
plot(time(Env1), signal(Env1),'r'); hold on 
plot(time(Env2), signal(Env2),'b'); 
plot(time(Env3), signal(Env3),'g');  Fs = Fs*10.31; [time1, signal1] = multiple_synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac);
plot(time1, signal1,'k'); ylabel('Amplitude'); xlabel('Time [s]'); title('Sub Nyquist Effect'); legend('Env 1', 'Env 2', 'Env 3', 'g(x)'); 
xlim([0,8e-3]); 

%--- Plotting PSD Graphs - different methods
figure(2); 
[pxx,w] = periodogram(signal, rectwin(length(signal)), length(signal), Fs,'PSD');
subplot(1,3,1); semilogy(w,rdivide(pxx,w)); title('PSD - FFT'); ylabel('volts ^2/Hz'); xlabel('Hz'); 
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time, signal,'mean','linear'); 
subplot(1,3,2); power = LSSA_helper(CScoeff,'PSD',Omega); 
semilogy(Omega,power); title('PSD - LSSA'); ylabel('volts ^2/Hz'); xlabel('Hz'); 
[acf, lags] = autocorr(signal, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(1,3,3); 
[freq, amp1] = my_FFT(time,acf); 
semilogy(freq,rdivide((amp1*2),1)); title('PSD - ACF'); ylabel('volts ^2/Hz'); xlabel('Hz'); 

%% Test 3 (one gap)
close all; clearvars; 
%--- Constructign signal with lin trend, var freq
amp = [10 20 420 1000 500]; 
phase = [0 60 20 15 10]; % degrees 
freq_start = [300 430 550 3000 4000]; freq_end = [300 430 750 3000 4000]; 
Fs = 10000; T = 2;  SNR_dB = 25;
linear_trend_frac = repmat(1,[1,5]); 
[time, signal] = multiple_synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac);

%--- Removing gap of len 450
[time1,signal1] = my_gappy(time,signal,8000,1000); 

%--- Interpolating (Linear, Quadractic, Cubic Spline, Sinc) 
Interp_lin = interp1(time1,signal1,time); 
Interp_spline = interp1(time1,signal1,time,'spline'); 
Interp_fft = interpft(signal1, length(time)); 

figure(1)
subplot(2,2,1); plot(time,Interp_lin); title('Lin Interp'); 
subplot(2,2,2); plot(time,signal); title('Original'); 
subplot(2,2,3); plot(time,Interp_fft); title('FFT Interp'); 
subplot(2,2,4) ; plot(time,Interp_spline); title('Spline Interp'); 

%--- Computing Power Spectrum 
figure(2); 
%--- Linear Interp
[pxx,w] = periodogram(Interp_lin, rectwin(length(signal)), length(signal), Fs,'power');
subplot(2,2,1); plot(w,10*log10(pxx)); title('Spec. - FFT - lin'); ylabel('dB'); xlabel('Hz')
[acf, lags] = autocorr(Interp_lin, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(2,2,2); 
[freq, amp] = my_FFT(time,acvf); 
plot(freq,10*log10((amp*2)/20)); title('Spec. - ACVF - lin'); ylabel('dB'); xlabel('Hz')

%---Spline Interp
[pxx,w] = periodogram(Interp_spline, rectwin(length(signal)), length(signal), Fs,'power');
subplot(2,2,3); plot(w,10*log10(pxx)); title('Spec. - FFT - spline'); ylabel('dB'); xlabel('Hz')
[acf, lags] = autocorr(Interp_spline, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(2,2,4); 
[freq, amp] = my_FFT(time,acvf); 
plot(freq,10*log10((amp*2)/20)); title('Spec. - ACVF - spline'); ylabel('dB'); xlabel('Hz')

%--- FFT Interp and LSSA Power
figure(3)
[pxx,w] = periodogram(Interp_fft, rectwin(length(signal)), length(signal), Fs,'power');
subplot(1,3,1); plot(w,10*log10(pxx)); title('Spec. - FFT - (FFT)'); ylabel('dB'); xlabel('Hz')
[acf, lags] = autocorr(Interp_lin, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(1,3,2); 
[freq, amp] = my_FFT(time,acvf); 
plot(freq,10*log10((amp*2))); title('Spec. - ACVF - (FFT)'); ylabel('dB'); xlabel('Hz')
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time, signal,[1, 1, 0, 0]); 
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time1, signal1,Omega,[1, 1, 0, 0]); 
subplot(1,3,3); power = LSSA_helper(CScoeff,'power'); plot(Omega,10*log10(power/4)); title('Spec - LSSA'); ylabel('dB'); xlabel('Hz')

%% Test 4 (two gap)
close all; clearvars; 
%--- Constructign signal with lin trend, var freq
amp = [10 20 420 1000 500]; 
phase = [0 60 20 15 10]; % degrees 
freq_start = [300 430 550 3000 4000]; freq_end = [300 430 750 3000 4000]; 
Fs = 10000; T = 2;  SNR_dB = 25;
linear_trend_frac = repmat(1,[1,5]); 
[time, signal] = multiple_synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac);

%--- Removing gap of len 450
[time1,signal1] = my_gappy(time,signal,8000,1000); 
[time1,signal1] = my_gappy(time1,signal1,14000,1000); 

%--- Interpolating (Linear, Quadractic, Cubic Spline, Sinc) 
Interp_lin = interp1(time1,signal1,time); 
Interp_spline = interp1(time1,signal1,time,'spline'); 
Interp_fft = interpft(signal1, length(time)); 

figure(1)
subplot(2,2,1); plot(time,Interp_lin); title('Lin Interp'); 
subplot(2,2,2); plot(time,signal); title('Original'); 
subplot(2,2,3); plot(time,Interp_fft); title('FFT Interp'); 
subplot(2,2,4) ; plot(time,Interp_spline); title('Spline Interp'); 

%--- Computing Power Spectrum 
figure(2); 
%--- Linear Interp
[pxx,w] = periodogram(Interp_lin, rectwin(length(signal)), length(signal), Fs,'power');
subplot(2,2,1); plot(w,10*log10(pxx)); title('Spec. - FFT - lin'); ylabel('dB'); xlabel('Hz')
[acf, lags] = autocorr(Interp_lin, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(2,2,2); 
[freq, amp] = my_FFT(time,acvf); 
plot(freq,10*log10((amp*2)/20)); title('Spec. - ACVF - lin'); ylabel('dB'); xlabel('Hz')

%---Spline Interp
[pxx,w] = periodogram(Interp_spline, rectwin(length(signal)), length(signal), Fs,'power');
subplot(2,2,3); plot(w,10*log10(pxx)); title('Spec. - FFT - spline'); ylabel('dB'); xlabel('Hz')
[acf, lags] = autocorr(Interp_spline, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(2,2,4); 
[freq, amp] = my_FFT(time,acvf); 
plot(freq,10*log10((amp*2)/20)); title('Spec. - ACVF - spline'); ylabel('dB'); xlabel('Hz')

%--- FFT Interp and LSSA Power
figure(3)
[pxx,w] = periodogram(Interp_fft, rectwin(length(signal)), length(signal), Fs,'power');
subplot(1,3,1); plot(w,10*log10(pxx)); title('Spec. FFT - (FFT)'); ylabel('dB'); xlabel('Hz')
[acf, lags] = autocorr(Interp_lin, 'NumLags', length(signal)-1); acvf = acf*rms(signal)^2; subplot(1,3,2); 
[freq, amp] = my_FFT(time,acvf); 
plot(freq,10*log10((amp*2))); title('Spec. ACVF - (FFT)'); ylabel('dB'); xlabel('Hz')
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time, signal,[1, 1, 0, 0]); 
[spectrum, CritVal, CScoeff, res, norm_res, coeff, cov, Omega] = LSSA(time1, signal1,Omega,[1, 1, 0, 0]); 
subplot(1,3,3); power = LSSA_helper(CScoeff,'power'); plot(Omega,10*log10(power/4)); title('Spec LSSA'); ylabel('dB'); xlabel('Hz')

