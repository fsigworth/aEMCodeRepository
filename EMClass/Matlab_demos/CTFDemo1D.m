% CTFDemo1D.m

clf;
slFontSize=18;
slx=.03;
sly=.85;
n=256;
g=randn(n,1);
fmax=0.3;
delta=.5;
f=(0:n-1)'*fmax/n;
x=(0:n-1)'/fmax/2;

subplot(3,2,1);
plot(x,[g 0*g]);
axis([0 max(x) -3 3]);
xlabel('x coordinate, Å');
ylabel('Signal');
SubplotLabel(slx,sly,'A',slFontSize);
title('Unfiltered Signal');

G=fftshift(fft(g));

subplot(3,2,2);
plot(f,abs(G).^2/n);
hold on;
plot(f,0*G+1,'k-');
hold off;
axis([0 fmax 0 5]);
xlabel('Frequency, Å^{-1}');
ylabel('Spectral density');
SubplotLabel(slx,sly,'B',slFontSize);

h=ContrastTransfer(f,.025,delta,2,10,.1);

subplot(3,2,4);
plot(f,[h nan*h 0*h]);
xlabel('Frequency, Å^{-1}');
ylabel('CTF');
SubplotLabel(slx,sly,'C',slFontSize);

Gc=G.*h;
gc=real(ifftn(ifftshift(Gc)));

subplot(3,2,5);
plot(x,gc);
axis([0 max(x) -1 1]);
xlabel('x coordinate, Å');
ylabel('Filtered Signal');
SubplotLabel(slx,sly,'D',slFontSize);
title('Filtered Signal');

subplot(3,2,6);
plot(f,abs(Gc).^2/n);
hold on;
plot(f,h.^2,'k-');
hold off;
axis([0 fmax 0 2.5]);
xlabel('Frequency, Å^{-1}');
ylabel('Spectral density');
SubplotLabel(slx,sly,'E',slFontSize);
