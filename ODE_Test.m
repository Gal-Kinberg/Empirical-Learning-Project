close all
clear
% clc

PLOT_MOVIE = false;
%%
L  = 5;
g  = 9.8; 
w0 = sqrt(g / L);
T  = 2 * pi / w0;
f0 = 1 / T;

Tmax = 35;
dt   = 0.1;
vT   = (0 : dt : Tmax)';
N    = length(vT);

%-- Initial Conditions: y0(1) is angle, y0(2) is angular velocity, y0(3) is
%-- pole length (L)
y0    = [pi/5 0 L];
damp  = 0.1;
ODE   = @(t,y) [y(2);
                -g / y(3) * sin(y(1)) - damp * y(2);
                0];
[~, mY] = ode45(ODE, vT, y0);

%%
if PLOT_MOVIE == true
    figure();
    for ii = 2 : numel(vT)
        plot(0, 0, '.black', 'MarkerSize', 20); %-- plot black dot at center
        
        %-- plot "mass":
        hold on; plot(mY(ii,3) * sin(mY(ii,1)), -mY(ii,3) * cos(mY(ii,1)), '.', 'MarkerSize', 45);
        set(gca, 'FontSize', 16);
        
        %-- plot "pole":
        plot([0 mY(ii,3) * sin(mY(ii,1))], [0 -mY(ii,3) * cos(mY(ii,1))], 'LineWidth', 3); hold off;
        
        %-- aesthetics:
        xlim([-1.5 * mY(ii,3) 1.5*mY(ii,3)]); ylim([-1.5*mY(ii,3) 1.5*mY(ii,3)]); grid on;
        title(['L = ',num2str(mY(ii,3)),'[m]  f_0 = ',num2str(f0),' [Hz]']);
        xlabel(['t = ',num2str(vT(ii)),' [sec]']);
        drawnow;
    end
end

%% Phone Data
% PhoneData = load('Merged_raw_data_3.csv');
% dt = 0.02;
% mY = PhoneData.


%% Diffusion Map
Fs = 1 / dt;
        
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
vYlim = ylim;
plot([f0, f0], [vYlim(1), vYlim(2)], ':r', 'LineWidth', 2 );

%% XY Diffusion

mY_XY = [L * sin(mY(:,1)), -L * cos(mY(:,1))];

mW         = squareform( pdist(mY_XY) );
eps        = median(mW(:));
mK         = exp(-mW.^2 / eps^2);
mA         = mK ./ sum(mK, 2);

N = size(mY_XY,1);
[mPhi, mLam] = eig(mA);
f            = Fs / 2 * linspace(-1, 1, N + 1); f(end) = [];

figure; hold on; set(gca, 'FontSize', 16);
plot(f, fftshift( abs( fft(mPhi(:,2)) ) ), 'LineWidth', 2 );
xlabel('f [Hz]'); title('Fourier of first (non-trivial) eigenvector, mY_{XY}');
vYlim = ylim;
plot([f0, f0], [vYlim(1), vYlim(2)], ':r', 'LineWidth', 2 );

%% Velocity Diffusion
vel = diff(mY_XY) / dt;
N = size(vel, 1);
mY_V = [ mY_XY(1:N,:) vel ];

mW         = squareform( pdist(mY_V) );
eps        = median(mW(:));
mK         = exp(-mW.^2 / eps^2);
mA         = mK ./ sum(mK, 2);


[mPhi, mLam] = eig(mA);
f            = Fs / 2 * linspace(-1, 1, N + 1); f(end) = [];

figure; hold on; set(gca, 'FontSize', 16);
plot(f, fftshift( abs( fft(mPhi(:,2)) ) ), 'LineWidth', 2 );
xlabel('f [Hz]'); title('Fourier of first (non-trivial) eigenvector, mY_V');
vYlim = ylim;
plot([f0, f0], [vYlim(1), vYlim(2)], ':r', 'LineWidth', 2 );

%% Plot with Phase Space Test
% figure();
% for i = 2:numel(t)
%     subplot(1,2,1);
%     plot(y(i,3)*sin(y(i,1)),-y(i,3)*cos(y(i,1)),'.','MarkerSize',45);
%     hold on; plot([0 y(i,3)*sin(y(i,1))],[0 -y(i,3)*cos(y(i,1))],'LineWidth',3); hold off;
%     xlim([-L L]); ylim([-1.5*L L]); grid on;
%     % xlabel t; ylabel y ;
%     title 'Real Space';
%     drawnow;
%     subplot(1,2,2);
%     hold on;
%     plot(y(i,1),y(i,2),'.b','MarkerSize',10);
%     hold on; plot([y(i-1,1) y(i,1)],[y(i-1,2) y(i,2)],'LineWidth',3); hold off;
%     xlim([min(y(:,1)) max(y(:,1))]); ylim([min(y(:,2)) max(y(:,2))]); grid on;
%     xlabel y; ylabel 'dy/dt';
%     title 'The Matrix';
%     drawnow;
