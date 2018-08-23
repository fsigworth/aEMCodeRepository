function sp=ccdAliasedGaussians(n,sigmas,amps)
% Compute the sum of 2D Gaussians having the standard deviations sigmas
% (normalized to total width n) and the amplitudes amps.  The value amps=1
% gives a Gaussian with the maximum value 1 in the center.  But, we assume
% unlimited aliasing (implemented by undersampling followed by FT), so
% sigmas larger than about .2 yield aliasing effects at the edges and
% increased value in the center.  Sigma >.7 yields essentially constant
% functions. This function is used in modeling camera noise spectra as sums of
% Gaussians.
sigmaXs=1./(2*pi*sigmas);  % real space sigmas
maxSx=max(sigmaXs);
if numel(n)<2
    n=n*[1 1];
end;

% We will compute real-space Gaussians on a region nx points square.
nx=n;  % Worst case: use the entire space.
if all(maxSx<n/12) % We're beyond 6 sigma (error < 2e-8) can truncate
    nx=[1 1]*ceil(maxSx*12);  % size of the truncated kernel
end;

hx=zeros(nx);  % real-space kernel
f2=double(Radius(nx)).^2;
for k=1:numel(sigmas)
    s=1./(2*pi*sigmas(k));
    a=2*pi*sigmas(k).^2;  % corresponding real-space sigma (not normalized)
    hx=hx+amps(k)*a*exp(-f2/(2*s^2));
end;
% Now pad the kernel up to the original size
if n>nx
    hx=Crop(hx,n,0,0);  % pad the kernel, force to be double
end;
sp=fftshift(real(fftn(ifftshift(hx))));


% % Old brute-force code that does the same thing for 1 Gaussian.  na is the
% % number of aliased copies of the function on each side.
% function sp=AliasedModelSpectrum3(n,f0,na)
% ct=ceil(n+1)/2;
% sp=zeros(n,n);
% for i=1:2*na+1
%     x=ct+(i-1-na)*n;
%     for j=1:2*na+1
%         y=ct+(j-1-na)*n;
%         f=RadiusNorm(n,[x,y]);
%         sp=sp+exp(-f.^2/(2*f0.^2));
%     end;
% end;
