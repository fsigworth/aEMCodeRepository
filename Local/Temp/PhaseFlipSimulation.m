% PhaseFlipSimulation.m

n=1024;
xs=(-n/2:n/2-1)'/32;
nr=6;
nc=1;

f=2;
sigma=.5;

img=cos(2*pi*f*xs).*exp(-xs.^2/2/sigma^2);
fimg=fftshift(fft(ifftshift(img)));
figure(1);
subplot(nr,nc,1);
plot(img);
ylabel('img');
subplot(nr,nc,2)
plot(real(fimg));
ylabel('fimg');
fc=2;
ctf=sin(2*pi*fc*xs);

fimg2=ctf.*fimg;

subplot(nr,nc,3);
plot(real(fimg2))
ylabel('fimg2');

subplot(nr,nc,4);
img2=imag(fftshift(ifft(ifftshift(fimg2))));
plot(img2);
ylabel('img2');

flip=sign(ctf);
% flip=ctf;
imgf=real(fftshift(ifft(ifftshift(fimg2.*flip))));
subplot(nr,nc,5);
plot(imgf)
% axis([-inf inf -0.5 0.5]);
ylabel('imgf');

% flip amplitude:
ampFlip=imgf(:)'*img(:)/sqrt(img(:)'*img(:))

% correl amplitude
ampNoFlip=img2(:)'*img2(:)/sqrt(img2(:)'*img2(:))

snrRatio=(ampFlip/ampNoFlip)^2
