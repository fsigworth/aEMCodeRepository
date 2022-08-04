function m=isft(fm);
% Shifted inverse, n-dimensional FFT
m=fftshift(ifftn(fftshift(fm)));
