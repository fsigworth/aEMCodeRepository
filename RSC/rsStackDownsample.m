function [siOut,imgsOut]=rsStackDownsample(si, imgs, nOut, mask)
% function [siOut,imgsOut]=rsStackDownsample(si, imgs, nOut, mask)
% Takes a StackInfo structure and stack of images and resamples them to
% the image size nOut.  The ctfs are correspondingly cropped as well.
% mask is an optional Fourier mask nOut in size, e.g.
%   nMsk=min(nIn,nOut);
%   mask=fuzzymask(nOut,2,nMsk*0.475,nMsk*.05);  % 90% of Nyquist
% which is used both in resampling the images and masking the ctfs.

if nargin<4
    mask=1;
end;
nImgs=size(imgs,3);
if nargout>1
    imgsOut=DownsampleGeneral(imgs,nOut,[],nImgs>1);
end;
ds=size(imgs,1)/nOut;
siOut=si;

siOut.pixA=si.pixA*ds;
siOut.yClick=si.yClick/ds;
siOut.rVesicle=si.rVesicle/ds;
if isfield(si,'mbnOffset')
    siOut.mbnOffset=si.mbnOffset/ds;
end;
siOut.ctfs=Crop(si.ctfs,nOut,1).*mask;
if isfield(si,'ctfGroups')
    siOut.ctfGroups=Crop(si.ctfGroups,nOut,1).*mask;
end;

% nc=size(si.ctfs,3);
% 
% for i=1:nc
%     siOut.ctfs(:,:,i)=Crop(si.ctfs(:,:,nc),nOut);
% end;
