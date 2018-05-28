%% Camera Test
close all
clear

dirPath = './Pandulum_Movies/';

%% Movies Description
%{
    27.5 Movies:
    Pend (1-4)      - Two Frequencies Tests
    Pend (5)        - 2D Frequncies Test
    Pend (6-7)      - No Phone Tests
    Pend (8)        - Simple Test
    Pend (9)        - Fail
    Pend (10)       - Dark Shot
    Pend (11-13)    - Close Angles
    Pend (14-15)    - Spring Tests
    Pend (16-17)    - Spring
    Pend (18-20)    - Pend + Spring
    Pend (21)       - Fluorescence
    Pend (22-23)    - Fluorescence Slo-Mo + Pend

    Pend_01         - Light Shot
    Pend_02-04      - Small, Mid, Large Amplitudes
    Pend_05         - Angle Shot
    Pend_06         - Angle Large Amp
    Pend_07         - Steady Hand-Held
    Pend_08         - Unstable Hand-Held
    Pend_09         - Tracking Hand-Held

    Pend_Self (1-4) - Shots From The Pendulum-Phone
%}
%% Read Movie

mov               = VideoReader('Pend (21).MP4');
mov.CurrentTime   = 0;             %-- Movie start time in seconds
EndTime           = mov.Duration;  %-- mov.Duration for whole movie
Fs                = mov.FrameRate;
dt                = 1 / Fs;
Scale             = 0.1;
video             = [];

%-- Cropped is for movies 21-24
while ( hasFrame(mov) && mov.CurrentTime <= EndTime )
    frame          = readFrame(mov);
%     cropped        = imcrop(frame,[900 290 500 450]);
%     rotated        = imrotate(cropped,90);
    scaled_frame   = imresize( rgb2gray(frame) ,Scale);
    video          = [video scaled_frame(:)];
end
mY = double(video');
%% Diffusion Map

mW         = squareform( pdist(mY) );
eps        = median(mW(:));
mK         = exp(-mW.^2 / eps^2);
mA         = mK ./ sum(mK, 2);

N = size(mY,1);
[mPhi, mLam] = eig(mA);
f            = Fs / 2 * linspace(-1, 1, N + 1); f(end) = [];

figure; hold on; set(gca, 'FontSize', 16);
plot(f, fftshift( abs( fft(mPhi(:,2)) ) ), 'LineWidth', 2 );
xlabel('f [Hz]'); title('Fourier of first (non-trivial) eigenvector, mY');
grid on;

% vYlim = ylim;
% plot([f0, f0], [vYlim(1), vYlim(2)], ':r', 'LineWidth', 2 );
