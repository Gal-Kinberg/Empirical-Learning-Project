close all
clear
% clc

%%
L=5; g=9.8; 
w0 = sqrt(g/L); T = 2*pi/w0; f0 =1/T;

Tmax = 35; dt = 0.1;
tspan = 0:dt:Tmax;

% Initial Conditions: y0(1) is angle, y0(2) is angular velocity, y0(3) is
% pole length (L)
y0 = [pi/5 0 L];
damp = 0.1
[t,y] = ode45(@(t,y) [(y(2)) ; -(g/y(3))*sin(y(1))-damp*y(2) ; 0], tspan, y0);


figure();
for i = 2:numel(t)
    plot(0,0,'.black','MarkerSize',20); % plot black dot at center
    % plot "mass":
    hold on; plot(y(i,3)*sin(y(i,1)),-y(i,3)*cos(y(i,1)),'.','MarkerSize',45);
    % plot "pole":
    plot([0 y(i,3)*sin(y(i,1))],[0 -y(i,3)*cos(y(i,1))],'LineWidth',3); hold off;

    %   aesthetics:
    xlim([-1.5*y(i,3) 1.5*y(i,3)]); ylim([-1.5*y(i,3) 1.5*y(i,3)]); grid on;
    title(['L = ',num2str(y(i,3)),'[m]  f0 = ',num2str(f0),' [Hz]']);
    xlabel(['t = ',num2str(t(i)),' [sec]']);
    drawnow;
end

%% Diffusion Map
Fs = 1/dt;

distances = squareform(pdist(y));
eps = median(y(:,1));
W = exp(-(distances.^2)./((100*eps)^2));
A = W./sum(W,2);
[Phi, Lambda] = eig(A);
f = linspace(-Fs/2,Fs/2,numel(y(:,1)));

figure;
stem(f,abs(fftshift(fft(Phi(:,2)))));
xlabel('f [Hz]'); title('Fourier of first (non-trivial) eigenvector');

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
% end