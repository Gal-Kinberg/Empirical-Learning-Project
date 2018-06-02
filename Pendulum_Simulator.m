close all
clear

PLOT_MOVIE = false;
%% Create Pendulum Simulation
L  = 5;
g  = 9.8; 
w0 = sqrt(g / L);
T  = 2 * pi / w0;
f0 = 1 / T;

Tmax = 35;
dt   = 0.1;
Fs   = 1 / dt;
vT   = (0 : dt : Tmax)';

%-- Initial Conditions: y0(1) is angle, y0(2) is angular velocity, y0(3) is
%-- pole length (L)
y0    = [pi/5 0 L];
damp  = 0.1;
ODE   = @(t,y) [y(2);
                -g / y(3) * sin(y(1)) - damp * y(2);
                0];
[~, mY] = ode45(ODE, vT, y0);

%% Plot Simulation Movie
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

%% Diffusion Map
        
[mPhi, mLam] = DiffusionMap(mY);

figure;
DiffusionPlot(mPhi,1,Fs,f0);

%% XY Diffusion

mY_XY = [L * sin(mY(:,1)), -L * cos(mY(:,1))];

[mPhi, mLam] = DiffusionMap(mY_XY);

figure;
DiffusionPlot(mPhi,1,Fs,f0);

%% Velocity Diffusion
vel = diff(mY_XY) / dt;
N = size(vel, 1);
mY_V = [ mY_XY(1:N,:) vel ];

[mPhi, mLam] = DiffusionMap(mY_V);

figure;
DiffusionPlot(mPhi,1,Fs,f0);
