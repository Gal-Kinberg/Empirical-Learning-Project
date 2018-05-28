close all
clear

PLOT_MOVIE = false;

dirPath = './Phone_Data/Phone_27.5/';

%% Movies Description
%{
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
    Pend (21)       - Fluorescent
    Pend (22-23)    - Fluorescent Slo-Mo + Pend

    Pend_01         - Light Shot
    Pend_02-04      - Small, Mid, Large Amplitudes
    Pend_05         - Angle Shot
    Pend_06         - Angle Large Amp
    Pend_07         - Steady Hand-Held
    Pend_08         - Unstable Hand-Held
    Pend_09         - Tracking Hand-Held

    Pend_Self (1-4) - Shots From The Pendulum Phone
%}

%% Load Phone Data
Fs        = 50;
dt        = 1 / Fs;
L         = 1.5;

fileName  = 'Pend_04.xlsx';
% [M, ~, X] = xlsread([dirPath, fileName]);
T         = readtable([dirPath, fileName]);

%% Extract Data
mY          = T.Acceleration_z_Cm_sec_2;
[~, minind] = min(mY);
mY          = mY(minind:end);
plot(mY);

%% Plot Simulation Movie
if PLOT_MOVIE == true
    figure();
    for ii = 2 : numel(mY)
        plot(0, 0, '.black', 'MarkerSize', 20); %-- plot black dot at center
        
        %-- plot "mass":
        hold on; plot(L * sin(mY(ii,1)), -L * cos(mY(ii,1)), '.', 'MarkerSize', 45);
        set(gca, 'FontSize', 16);
        
        %-- plot "pole":
        plot([0 L * sin(mY(ii,1))], [0 -L * cos(mY(ii,1))], 'LineWidth', 3); hold off;
        
        %-- aesthetics:
        xlim([-1.5 * L 1.5*L]); ylim([-1.5*L 1.5*L]); grid on;
%         title(['L = ',num2str(mY(ii,3)),'[m]  f_0 = ',num2str(f0),' [Hz]']);
%         xlabel(['t = ',num2str(vT(ii)),' [sec]']);
        drawnow;
    end
end

%% Diffusion Map
        
mW         = squareform( pdist(mY) );
eps        = median(mW(:));
mK         = exp(-mW.^2 / eps^2);
mA         = mK ./ sum(mK, 2);

% N = size(mY,1);
N = 2^12;
[mPhi, mLam] = eig(mA);
f            = Fs / 2 * linspace(-1, 1, N + 1); f(end) = [];

figure; hold on; set(gca, 'FontSize', 16);
plot(f, fftshift( abs( fft( ( mPhi(:,2) ), N ) ) ), 'LineWidth', 2 );
xlabel('f [Hz]'); title('Fourier of First (non-trivial) Eigenvector');
grid on;
% vYlim = ylim;
% plot([f0, f0], [vYlim(1), vYlim(2)], ':r', 'LineWidth', 2 );
