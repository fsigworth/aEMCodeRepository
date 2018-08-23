function [projs angles]=spReproject(vol,angStep)
% function [projs angles]=spReproject(vol,angStep)
% Create a set of projections of the volume vol, and return the [phi theta]
% angles corresponding to the projection directions.
ovs=4; % oversampling factor

% make the projection angles
npsi=ceil(360/angStep);
ntheta=round(npsi/4);
angs=SphereAngles2(ntheta,npsi);
angles=angs*180/pi;  % in degrees

% Make the oversampled Fourier volume
nx=size(vol,1)*ovs;
rxvol=Crop(vol,nx);  % real-space padded volume
fxvol=fftshift(fftn(ifftshift(rxvol))); % oversampled Fourier volume

slices=OversamplingProject2(fxvol,angs,ovs);
projs=fftshift(real(ifft2(ifftshift(slices))));
