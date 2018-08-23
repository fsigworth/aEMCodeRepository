function sp=AliasedGaussians(n,sigmas,amps)
% Compute the sum of 2D Gaussians having the standard deviations sigmas
% (normalized to total width n) and the amplitudes amps.  The value amps=1
% gives a Gaussian with nominal value 1 in the center.  But, we assume
% unlimited aliasing, so sigmas larger than about .2 yield aliasing effects
% at the edges.  This function is used in modeling camera spectra.
h=zeros(n,n);
f2=double(Radius(n)).^2;
        for k=1:numel(sigmas)
            s=1/(2*pi*sigmas);
            a=2*pi*sigmas^2;
            h=h+amps(k)*a*exp(-f2/(2*s^2));
        end;
        sp=fftshift(real(fftn(ifftshift(h))));
   sp3=sp;     
imags(sp)
% plot(sp)