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

mov               = VideoReader('Pend (11).MP4');
mov.CurrentTime   = 0;             %-- Movie start time in seconds
EndTime           = mov.Duration;  %-- mov.Duration for whole movie
Fs                = mov.FrameRate;
dt                = 1 / Fs;
FramsNum          = Fs * mov.Duration;
Scale             = 0.1;
video             = [];

Chain = 10;
%-- Cropped is for movies 21-24
while ( hasFrame(mov) && mov.CurrentTime <= (EndTime - (Chain + 1) * dt) )
    frame          = readFrame(mov);
    %     cropped        = imcrop(frame,[900 290 500 450]);
    %     rotated        = imrotate(cropped,90);
    scaled_frame   = imresize( rgb2gray(frame) ,Scale);
    temp = scaled_frame;
%   chain several frames together:
    if Chain > 1
        for i = 1 : (Chain - 1)
            if(hasFrame(mov))
                frame          = readFrame(mov);
                scaled_frame   = imresize( rgb2gray(frame) ,Scale);
                temp           = [temp(:) ; scaled_frame(:)];
            end
        end
    end
%   add current frames to final matrix
    video              = [video temp];
end 
mY = double(video');
%% Diffusion Map

[mPhi, mLam] = DiffusionMap(mY);

figure;
DiffusionPlot(mPhi,1,Fs / Chain,0);
