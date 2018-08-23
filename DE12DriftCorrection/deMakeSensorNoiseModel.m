function [s f2]=deMakeSensorNoiseModel(s, doTest)
% function s=deMakeSensorNoiseModel(s, doTest)
% Given the data structure s for a stack of DE images, first find the
% obvious hot pixel outliers from the mean, and store them as
% s.hotPixelMap.  Then use the power beyond 1/2 Nyquist to estimate the
% remaining hot pixel noise spectrum and fit it.  The goal is to eliminate
% the false peak in cross-correlations. The model of the spurious
% cross-spectrum is a constant plus 2 Gaussians in the x direction, and one
% in the y direction. Two of these are very narrow (about 1 pixel).  The
% parameters are stored in s.modelPars. The model is
%  S(fx,fy)=amp1 + \sum _ j ampj * exp(-(f*p(j))^2/2)
% where f is the normalized frequency in x or y, according to the value of
% s.modelUseX(j).
% The fields stored in s are
%  s.hotPixelMap
%  s.modelPars
%  s.modelUseX
% The idea is that linear least-squares will be used to determine the
% amplitudes in each case.
% 
% The model assumes only a constant spectral density along the horizontal
% center line.  This is not a good model for that noise, and it's probably
% best to zero out that line in masking spectra.
%
%   fs April 2012

if nargin<2
    doTest=0;
end;
% 
% doTest=1;

displayOn=1;
dbin=4;  % display binning
nsds=8;  % threshold, in sds, for hot pixel mask
mskWidth=0.1;    % Fraction of Fourier space to mask out for variance estimate.
ds=16;  % downsampling for spectrum fitting

[nx ny nim]=size(s.imgs);
sz=[nx ny];
% Get the outlier map and store it as s.hotPixelMap
[im1 im2]=deGetAlternateMeans(s);
s.hotPixelMap=deHotPixelMap(im1+im2,nsds);

[msk mskx msky]=deMakeSpectrumMask(sz,mskWidth,ds);

%% Take means of alternate frames as representative of the
% noise.  Take its cross-spectrum with the rest of the movie frames.
im1=double(RemoveOutliers(im1,s.hotPixelMap));
im1x=im1-mean(im1(:));
im2=double(RemoveOutliers(im2,s.hotPixelMap));
im2x=im2-mean(im2(:));  % keep these 'x' images around for later.

%% Compute the cross-spectrum, and fit the spectrum of the hf noise.
crossSpec=fftshift(real(fftn(im1x).*conj(fftn(im2x))))/numel(im1x);

%% Fit the spectrum of the hf noise, which we take to be sensor noise.
% Fitting the full 2D spectrum takes too long, so we fit msk times the
% downsampled spectrum, and very selected regions (m0) time the
% full-resolution spectrum.

% operate on downsampled points
mskd=BinImage(msk,ds);
mskdx=BinImage(mskx,ds);  %% illustration
% Zero out the horizontal and vertical line
[nxd nyd]=size(mskd);
mskd(nxd/2+1,:)=0;  % vertical line
mskd(:,nyd/2+1)=0;  % horizontal line
mskdx(:,nyd/2+1)=0;  % horizontal line

csd=BinImage(crossSpec,ds);
m1=mskd(:);  % masked points
q1=find(m1);
c1=csd(q1); % data as vector
nq1=numel(q1);

% operate on original points
m0=zeros(nx,ny);
dw=4;  % pick up vertical and horizontal bands in the middle.
m0(nx/2+1-dw:nx/2+1+dw,:)=1;
m0(:,ny/2+1-dw:ny/2+1+dw)=1;
m0=m0.*msk;  % restrict the bands to the msk area
q0=find(m0);
c0=crossSpec(q0);

nxd=nx/ds;
nyd=ny/ds;
[xd yd]=ndgrid((-nxd/2:nxd/2-1)'/nxd,(-nyd/2:nyd/2-1)'/nyd);  % normalized freq coords
[x y]=ndgrid((-nx/2:nx/2-1)'/nx,(-ny/2:ny/2-1)'/ny);

xs=[xd(q1);x(q0)];
ys=[yd(q1);y(q0)];
cs=[c1;c0];  % data vector
weights=[0*c1+1;0*c0+.1];  % smaller weight on ufiltered data.

pvals=[2 100 4000]; % Gaussian decay parameter, constant
useX=[1 1 0];
% xvals=[xs xs ys];
np=numel(pvals);
p=Simplex('init',pvals);
fmat=zeros(numel(cs),np+1);
fmat(:,1)=weights;
model=zeros(nxd,nyd);
model0=model;

for i=1:200  % iterations
    for j=1:np
        if useX(j)
            fmat(:,j+1)=weights.*exp(-abs(p(j)*xs).^2/2);
        else
            fmat(:,j+1)=weights.*exp(-abs(p(j)*ys).^2/2);
        end;
    end;
    amps=LinLeastSquares(fmat,weights.*cs);
    d=fmat*amps-weights.*cs;  % error
        cdat=csd.*mskd;
    p=Simplex(d'*d);
    if mod(i,10)==0 && displayOn
        model0(q1)=fmat(1:nq1,:)*amps;
        model=model0.*mskd;
        subplot(221);
        plot([(sum(cdat)./sum(mskd))' (sum(model)./sum(mskd))']);
        subplot(222);
        plot([sum(cdat,2)./sum(mskd,2) sum(model,2)./sum(mskd,2)]);
        subplot(223);
        imacs(model-cdat);
        subplot(224);
        imacs(abs(model).^.2);
        title(num2str(([i p amps'])));
        drawnow;
    end;
end;
%%  Construct diagnostic variables and plot them.
        subplot(221);
        f2.xvals=(sum(csd.*mskdx)./sum(mskdx))';
        f2.xmodel=(sum(model0.*mskdx)./sum(mskdx))';
        plot([f2.xvals f2.xmodel]);
        axis([0 inf 0 1.1*max(f2.xvals)]);
        
        
        subplot(222);
        f2.yvals=sum(cdat,2)./sum(mskd,2);
        f2.ymodel=sum(model,2)./sum(mskd,2);
        plot([f2.yvals f2.ymodel]);
        subplot(223);
        
        f2.mskd=mskd;
        f2.csd=csd;
%%

s.modelPars=[0 p];  % include the constant component
s.modelUseX=[1 useX];
s.modelAmps=amps;

if doTest
[imgs1 imgs2]=deMakeCorrImageSet(s,0,0);
    
    %%  Test subtraction in the whole image using the model
    ind1=2;
    ind2=3;
    
    % compute a new cross-spectrum
% [im1 im2]=deGetCorrImagePair(s,ind,1);

im1=imgs1(:,:,ind1);
im2=imgs2(:,:,ind2);
im1x=double(im1-mean(im1(:)));
im2x=double(im2-mean(im2(:)));

% Compute the cross-spectrum, and fit the spectrum of the hf noise.
crossSpec=fftshift(real(fftn(im1x).*conj(fftn(im2x))))/numel(im1x);

    p=s.modelPars;
    np=numel(s.modelPars);
    
    bfmat=zeros(nx*ny,np);
    bmsk=msk(:);
 
    [x y]=ndgrid((-nx/2:nx/2-1)'/nx,(-ny/2:ny/2-1)'/ny);
    
    for j=1:np
        if j<np fvals=x(:); else fvals=y(:); end;
        bfmat(:,j)=exp(-abs(p(j)*fvals).^2/2);
    end;
    mbfmat=bfmat.*repmat(bmsk,1,np);
    bamps=LinLeastSquares(mbfmat,bmsk.*crossSpec(:));
    
    bigModel=reshape(bfmat*bamps,nx,ny);
%     subplot(221);
%     imacs(BinImage(msk.*bigModel-msk.*crossSpec,16));


    CorrCrossSpec=crossSpec-bigModel;
    cc=real(fftshift(ifftn(ifftshift(CorrCrossSpec))));
    subplot(221);
    imacs(Crop(cc,32));
    subplot(223);
    plot(sect(Crop(cc,32)));
    
    cc0=real(fftshift(ifftn(ifftshift(crossSpec))));
    subplot(222);
    imacs(Crop(cc0,32));
    subplot(224);
    plot(sect(Crop(cc0,32)));
     
    
end;
