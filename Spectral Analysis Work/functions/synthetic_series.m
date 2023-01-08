%% Summary of synthetic_series(*args)

% This function constructs a synthetic series for one sine wave
% Inputs are positional arguments as follows: 
% (1) Set ampitude value
% (2) Set frequency value (start)
% (3) Set frequency value (end). If different then start frequency, freq
%     will vary linearly with time to the end frequency. 
% (4) Set phase value (degrees)
% (5) Set SNR - interpreted as SNR = (Amplitude_Signal / Amplitude_Noise)^2
%             - thereby settign SNR sets noise level sinusoid expressed in
%             dB
% (6) Add linear trend fraction - interpreted as fraction of specified amplitude
% (7) duration of signal (T)
% (8) sample rate of series (L)

function [t, signal] = synthetic_series(amp, freq_start, freq_end, phase, T, Fs, SNR_dB, linear_trend_frac)
    %% Constructing Signal
    Delta = 1/Fs;    %Sample Period
    N = T/Delta;     %Duration L             
    t = linspace(0,T-Delta,N);   %Time Vector
    freq_vector = linspace(freq_start,freq_end,N); %constructing varying frequency vector
    signal_pure = amp*sin(2*pi*freq_vector.*t + deg2rad(phase));
    
    %% Constructing Noise (normal dist.) with mean 0 and std = sqrt(1)
    rng(0,'twister');
    noise = normrnd(0,sqrt(1),[1,length(signal_pure)]); 
    
    %% Scaling noise and corrupting signal based on desired SNR
    P_signal = rms(signal_pure)^2; 
    P_noise = rms(noise)^2;
    SNR = 10^(SNR_dB/10);
    scale = sqrt(P_signal/(P_noise*SNR));     
    signal = signal_pure + scale*noise; 
    
    %% adding a linear trend as a function of 25% amplitude of wave
    lin_trend = linspace(0,linear_trend_frac*amp,N);
    signal = signal + lin_trend;
    
    
end