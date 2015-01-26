clear all; close all

%Load data file
d = load('corrupt.dat');
w = d(:,1);
x = d(:,2);

%Gaussian (instrument response), centered on t=0
n = length(x);
g(1:n/2) = exp(-(1:n/2).^2/2048^2);
gt = [g fliplr(g)];
gt = gt./sum(gt);

%FFT data and instrument response
xfft = fft(hamming(n).*x);
figure(10); loglog(xfft(1:n/2).*conj(xfft(1:n/2))); 
gtfft = fft(gt.');
%Scale gtfft and plot on top of xfft
gtfft_sc = gtfft.*(max(abs(xfft))/max(abs(gtfft)));
figure(10); hold on; loglog(gtfft_sc(1:n/2).*conj(gtfft_sc(1:n/2)), 'r');
title('FFT of data and instrument response', 'FontSize', 18)
set(gca, 'FontSize', 14)
xlabel('Frequency (dB)', 'FontSize', 16)
ylabel('Power (dBm)', 'FontSize', 16)

%Deconvolve instrument response from input signal
%Try a range of cutoff frequencies
cut = 19;
cutrange = cut-4:cut+2
ncut = length(cutrange)
frange = 1:length(cutrange)
fidx = 1;
figure; subplot(ncut, 1, 1)
for cutoff = cutrange
    idx = [1:cutoff (n-cutoff+2):n];
    dc_fft = zeros(1,n);
    dc_fft(idx) = xfft(idx)./gtfft(idx);
    xr(fidx,:) = abs(ifft(dc_fft));
    subplot(ncut, 1, fidx); 
    plot(xr(fidx,:)); title(['Reconstructed signal ' num2str(cutoff)]);
    fidx = fidx + 1;
end

%Plot the change between reconstructions using a mean-squared metric
mse = mean((diff(xr).^2)');
figure; plot(cutrange(2:end), log10(mse), 'x', 'MarkerSize', 10)
set(gca, 'FontSize', 14)
title('Change in reconstructions (MSE)', 'FontSize', 18)
xlabel('Cutoff bin', 'FontSize', 16)
ylabel('Mean-squared difference from previous reconstruction', 'FontSize', 16)
axx = axis; axx([1 3]) = axx([1 3]) -1; axx([2 4]) = axx([2 4]) + 1; axis(axx)



