function [imgs1 imgs2 spDens]=deMakeCorrImageSet(s,doPW,doNoise)
% function [imgs1 imgs2]=deMakeCorrImageSet(s,doPW,doNoise) Given the data
% structure s for a stack of DE images, remove outliers and use the
% variance maps to compute
% anti-noise and add it to independently-corrected pairs of images
% imgs1 and imgs2.  By default, the images are then pre-whitened, according
% to the detector's white noise spectrum, and anti-noise is added.
% spDens is the average
% high-frequency spectral density for each of the imgs1 stack.
%
if nargin<2
    doPW=1;
end;
if nargin<3
    doNoise=1;
end;

% Make the pre-whitening filter
if doPW
    pwFilt=ifftshift(1./sqrt(CCDModelSpectrum2D(3))');  % full size DE camera spectrum
else
    pwFilt=1;
end;

pwFilt=double(pwFilt);

[nx ny nim]=size(s.imgs);
sz=[nx ny];
imgs1=single(zeros(nx,ny,nim));
imgs2=single(zeros(nx,ny,nim));

if doNoise
    sd1=sqrt(mean(s.varMap1(:)));
    sd2=sqrt(mean(s.varMap2(:)));
end;

if nargout>2
    % peripheral mask for estimating spectral density
    spectMask=Cropo(single(zeros(sz/2)),sz,0,1); % fill outer 3/4 of space with ones.
    spectMask=spectMask/(numel(spectMask)*sum(spectMask(:)));
    spDens=zeros(nim,1);
end;

for i=1:nim
    [im1 im2]=deGetCorrImagePair(s,i);
    im1=RemoveOutliers(im1,s.hotPixelMap);
    im2=RemoveOutliers(im2,s.hotPixelMap);
    
    if (~doNoise) && (~doPW)  % simply copy the images
        imgs1(:,:,i)=im1;
        imgs2(:,:,i)=im2;
        
    else  % We're doing noise or prewhitening
        
        if doNoise
            v1=Downsample(s.varMap1(:,:,i),sz);  % interpolate to full size
            v2=Downsample(s.varMap2(:,:,i),sz);
            % We do the large FTs in double precision
            im1d=double(im1+v1.*s.noise1/sd1+v2.*s.noise2/sd2);
            im2d=double(im2-sd1*s.noise1-sd2*s.noise2);
        else  % we'll do pw only, no noise added.
            im1d=double(im1);
            im2d=double(im2);
        end;
        
        fim1=fftn(im1d-mean(im1d(:))).*pwFilt;  % pwFilt=1 if no prewhitening.
        fim2=fftn(im2d-mean(im2d(:))).*pwFilt;
        if nargout>2
            spDens(i)=abs(fim1(:)').^2*spectMask(:);
        end;
        imgs1(:,:,i)=real(single(ifftn(fim1)));
        imgs2(:,:,i)=real(single(ifftn(fim2)));
    end
    
    % %   Verify that added anti-Noise is small.
    %   sp1=RadialPowerSpectrum(im1);
    %   sp2=RadialPowerSpectrum(im1d);
    %   sp3=RadialPowerSpectrum(imgs1(:,:,i));
    %   semilogy([sp1 sp2 sp3]);
    %   title(i);
    %   drawnow;
end;

return

%% test code

% look at cc of whole image
[imgs1 imgs2]=deMakeCorrImageSet(s,0,1);
cc=real(ifftn(fftn(imgs1(:,:,1)).*conj(fftn(sum(imgs2,3)-imgs2(:,:,1)))));
cc1=Crop(fftshift(cc),128);
subplot(223);
imacs(cc1);
[imgs10 imgs20]=deMakeCorrImageSet(s,0,0);
cc0=real(ifftn(fftn(imgs10(:,:,1)).*conj(fftn(sum(imgs20,3)-imgs20(:,:,1)))));
cc10=Crop(fftshift(cc0),128);
subplot(224);
imacs(cc10);
subplot(221);
plot(sect(cc0));
plot(sect(cc1));
subplot(222);
plot(sect(cc10));




% Look at cc of part of an image.
ind=5;
ncrp=4096*3/4;
ncrp=512;
i1=imgs1(:,:,ind);
i2=mean(imgs2(:,:,6:9),3);
% [i1 ix]=deGetCorrImagePair(s,ind);
% [ix i2]=deGetCorrImagePair(s,ind+1);
%       i2=(sum(imgs2,3)-imgs2(:,:,ind))/(nim-1);
%      [i1 i2]=deGetCorrImagePair(s,ind,0);
i1=Crop(i1,ncrp);
i2=Crop(i2,ncrp);
i1=i1-mean(i1(:));
i2=i2-mean(i2(:));
cc=fftshift(real(ifftn(fftn(i1).*conj(fftn(i2)))));
ccf=GaussFilt(Crop(cc,512),.1);
subplot(221);
imacs(BinImage(i1,8));
subplot(223);
imacs(Crop(ccf,128));
subplot(224);
plot(sect(Crop(ccf,128)),'.-','markersize',8);
drawnow;