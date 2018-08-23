function [vol normvol]=spReconstruct(projs,weights,angles,k)
[n ny nim]=size(projs);

angrs=angles*pi/180;  % to radians
slices=fftshift(fft2(ifftshift(projs)));

% slice0=fuzzymask(n,2,n/2-3,2);
% slices0=repmat(slice0,[1 1 nim]);
slices0=zeros(n,n,nim);
for i=1:nim
    slices0(:,:,i)=weights(i);
end;
normvol=OversamplingReconstr2(slices0,angrs);
dvol1=OversamplingReconstr2(slices,angrs);

% Wiener filter "normalization"
epsi=nim*k;

dvol=dvol1.*normvol./(epsi+normvol.^2);
vol=fftshift(real(ifftn(ifftshift(dvol))));
