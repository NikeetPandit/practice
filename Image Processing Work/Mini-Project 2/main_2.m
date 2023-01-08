%Assignment 2 by Nikeet Pandit
%GS/MATH 6920 Harmonic Analysis and Image Processing

% --- Question List --- %
% 4.9a) Write function for band-reject filtering in FREQ domain % SEE GitHub
%
% 4.9b) Write function band-pass in FREQ domain % SEE GitHub
%
% 4.9c) Generate ideal (square) band-reject filter 
%
% 4.9d) Generate gaussian band-reject filter
%
% 4.9e) Generate butterworth band-reject filter
%
% 4.9f) Repeat (c)-(e) w/ band-pass
%
% --- Question List -- % 

%% 4.9 (c) (d) (e)
clearvars; close all
%Plotting Responses (c) (d) (e) Band-Reject Filters
fig1 = figure(1); 
subplot(1,3,1); 
H = brFilterTF4e('ideal', 512, 512, 128,60, 0); %See FUN on GitHub
h = surf(H); set(h,'LineStyle','none'); %So Surf Plot has colors since mesh is too fine
title('Ideal BR'); axis tight; pbaspect([2 2 1]) %Controlling aspect ratio of 3D plot

subplot(1,3,2); 
H = brFilterTF4e('gaussian', 512, 512, 128,60, 0);
h = surf(H); title('Gaussian BR'); axis tight; pbaspect([2 2 1])
set(h,'LineStyle','none'); 
subplot(1,3,3); 
H = brFilterTF4e('butterworth', 512, 512, 128,60, 1);
h = surf(H); title('Butterworth BR'); 
set(h,'LineStyle','none'); 
set(fig1, 'Position',[1, 1, 1500, 600]); axis tight; 
colorbar; pbaspect([2 2 1])

%% 4.9 f)
%Plotting Responses (f) for (c) (d) (e) Band-Pass Fitlers 
fig2 = figure(2); 
subplot(1,3,1); 
H = bpFilterTF4e('ideal', 512, 512, 128,60, 0);
h = surf(H); set(h,'LineStyle','none'); title('Ideal BP'); axis tight; pbaspect([2 2 1])

subplot(1,3,2); 
H = bpFilterTF4e('gaussian', 512, 512, 128,60, 0);
h = surf(H); title('Gaussian BP'); axis tight; 
set(h,'LineStyle','none'); pbaspect([2 2 1])

subplot(1,3,3); 
H = bpFilterTF4e('butterworth', 512, 512, 128,60, 1);
h = surf(H); title('Butterworth BP'); axis tight;
set(h,'LineStyle','none'); 
set(fig2, 'Position',[1, 1, 1500, 600]);
colorbar; pbaspect([2 2 1])


