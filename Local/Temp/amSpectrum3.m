function sp=amSpectrum3(n,fs,amps)
% n=256;
% fs=.4;
% amps=1;

ct=ceil((n+1)/2);
h=zeros(n,n);
f2=double(Radius(n)).^2;
        for k=1:numel(fs)
            s=1/(2*pi*fs);
            a=2*pi*fs^2;
            h=h+amps(k)*a*exp(-f2/(2*s^2));
        end;
        sp=fftshift(real(fftn(ifftshift(h))));
   sp3=sp;     
imags(sp)
% plot(sp)