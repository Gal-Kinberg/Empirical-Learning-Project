function  DiffusionPlot(mPhi, VecNum, Fs, f0)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

N = 2^12;
f = Fs / 2 * linspace(-1, 1, N + 1); f(end) = [];

hold on; set(gca, 'FontSize', 16);
plot(f, fftshift( abs( fft(mPhi(:,VecNum + 1),N) ) ), 'LineWidth', 2 );
xlabel('f [Hz]'); title(['Fourier of Eigenvector Number ',num2str(VecNum)]);
grid minor;

if f0 ~= 0
    vYlim = ylim;
    plot([f0, f0], [vYlim(1), vYlim(2)], ':r', 'LineWidth', 2 );
end
end

