%% Video Simulator
close all
clear

SPRING = true;
PLOT_MOVIE = true;
RECORD_MOVIE = true;
%% Create Pendulum Simulation
L  = 5;
g  = 9.8;
wp = sqrt(g / L);
T  = 2 * pi / wp;
f0 = 1 / T;

Tmax = 35;
dt   = 0.1;
Fs   = 1 / dt;
vT   = (0 : dt : Tmax)';

%-- Initial Conditions: y0(1) is angle, y0(2) is angular velocity, y0(3) is
%-- pole length (L)
y0    = [pi/5 0 L];
damp  = 0.0;
ODE   = @(t,y) [y(2);
                -g / y(3) * sin(y(1)) - damp * y(2);
                0];
[~, mY] = ode45(ODE, vT, y0);

%% Create Spring Simulation
if SPRING == true
K    = 10;
m    = 1;
w_sp = sqrt(K / m);
T_sp = 2 * pi / w_sp;
f_sp = 1 / T_sp;

%-- Initial Conditions: y0(1) is hight, y0(2) is velocity
y0_sp    = [0.3*L 0];
damp_sp  = 0.0;
ODE_sp   = @(t,y) [y(2);
                   -w_sp^2 * y(1) - damp_sp * y(2)];
[~, mY_sp] = ode45(ODE_sp, vT, y0_sp);

%-- spring parameters:
xa = 0.8 * L;  ya = 0;  %-- base point
ne = 7;                 %-- number of coils
a  = L;                 %-- natural length
ro = 0.1;               %-- radius

[xs,ys] = spring(xa,ya,0,0,ne,a,ro); %-- spring initialization
end
%% Simulate Movie
F(1) = struct('cdata',[],'colormap',[]);
if PLOT_MOVIE == true
    figure();
    for ii = 1 : numel(vT)
        
        %-- plot black dot at center:
        plot(0, 0, '.black', 'MarkerSize', 20); 
        
        %-- plot "mass":
        hold on; plot(mY(ii,3) * sin(mY(ii,1)), -mY(ii,3) * cos(mY(ii,1)), '.', 'MarkerSize', 100);
        set(gca, 'FontSize', 16);
        
        %-- plot "pole":
        plot([0 mY(ii,3) * sin(mY(ii,1))], [0 -mY(ii,3) * cos(mY(ii,1))], 'LineWidth', 3); 
        
        if SPRING == true
        %-- plot spring mass
        plot( xa , -0.7 * L + mY_sp(ii,1), '.red', 'MarkerSize', 100);
        
        %-- plot spring
        plot(xa, ya, '.black', 'MarkerSize', 20); 
        [xs,ys] = spring(xa, ya, xa , -0.7 * L + mY_sp(ii,1)); 
        plot(xs,ys,'LineWidth',3);
        end
        
        hold off;
        
        %-- aesthetics:
        xlim([-L, L]);
        ylim([-1.5 * L, 0.1]); grid on;
%         title(['L = ',num2str(mY(ii,3)),'[m]  f_0 = ',num2str(f0),' [Hz]']);
%         xlabel(['t = ',num2str(vT(ii)),' [sec]']);
        xticks([]);
        yticks([]);
%         drawnow;
        box off;
        axis off;
        
        F(ii) = getframe(gcf);
    end
end
close;

if RECORD_MOVIE == true
v = VideoWriter('SimulationVideo.avi');
open(v); writeVideo(v, F); close(v);
end
%% Read Movie

mov   = VideoReader('SimulationVideo.avi');
Scale = 0.1;
video = [];
while hasFrame(mov)
    frame = imresize( rgb2gray( readFrame(mov) ), Scale);
    video = [video frame(:)];
end
mY = double(video');
%% Diffusion Map

[mPhi, mLam] = DiffusionMap(mY);

N = 2^12;
f = Fs / 2 * linspace(-1, 1, N + 1); f(end) = [];

figure;
if SPRING == true
    subplot(121);
end
hold on; set(gca, 'FontSize', 16);
plot(f, fftshift( abs( fft(mPhi(:,2), N ) ) ), 'LineWidth', 2 );
xlabel('f [Hz]'); title('Fourier of First (non-trivial) Eigenvector');
vYlim = ylim;
plot([f0, f0], [vYlim(1), vYlim(2)], ':r', 'LineWidth', 2 );
grid minor;

if SPRING == true
    plot([f_sp, f_sp], [vYlim(1), vYlim(2)], ':g', 'LineWidth', 2 );
    legend('Fourier Transform','f_{pend}','f_{spring}');
    
    subplot(122); hold on; set(gca, 'FontSize', 16);
    plot(f, fftshift( abs( fft(mPhi(:,3), N) ) ), 'LineWidth', 2 );
    xlabel('f [Hz]'); title('Fourier of Second (non-trivial) Eigenvector');
    vYlim = ylim;
    plot([f0, f0], [vYlim(1), vYlim(2)], ':r', 'LineWidth', 2 );
    grid minor;
    plot([f_sp, f_sp], [vYlim(1), vYlim(2)], ':g', 'LineWidth', 2 );
    legend('Fourier Transform','f_{pend}','f_{spring}');
end