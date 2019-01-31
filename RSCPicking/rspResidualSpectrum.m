function [sval,xvals] = rspResidualSpectrum(mi,rscc,dis,coords,displayOn)
% Returns the maximum of the power spectrum of the vicinity of the
% particle.  If called with displayOn, makes a display of the residuals in
% Figure 2.
% The output argument xvals is a row vector of absolute values of extrema
% after various Gaussian filters.
%
% The spectrum is computed from the image as processed by the pre-whitening
% filter used in rsPickingPreprocessor to create model particles. The 
% weighting of spectral densities is by a Gaussian with sigma hard-wired
% to 2 (!) with the frequency step set by nd, the model particle box size.
xFreqs=[.02 .02 .01 .01]; % in inverse Å, equiv to sds of 4, 8, 16, 32QQQ Å

scaleUp=1;
sval=0;
pixA=mi.pixA*dis.ds;

if ~isfield(rscc,'eigenImgs') || size(rscc.eigenImgs,1)<2
    return
end;

partMaskRadiusA=dis.spectrumMaskRadiusA; % mask of particle for residual calculation
partMaskRadius=partMaskRadiusA/(dis.ds*mi.pixA); % radius in pixels
surrMaskRadius=partMaskRadius*1.4; % Surround mask for DC calculation
% cropBoxSize=NextNiceNumber(partMaskRadius*2.6);

coords=coords(:)';
ix=round(coords(1)/dis.ds+1);
iy=round(coords(2)/dis.ds+1);
pos=[ix iy];
% sval=0;
gfc=.14;  % cutoff of Gauss filter for display
nd=size(rscc.eigenImgs,1);
particleMask=fuzzymask(nd,2,partMaskRadius,partMaskRadius/2);
surroundMask=fuzzymask(nd,2,surrMaskRadius,partMaskRadius/3)-particleMask; % annular mask

% Get the theoretical particle image, undoing the pre-whitening that the
% picking preprocessor did.
filtParticle=scaleUp*rspMakeModelParticle(dis,mi,coords,rscc);
% Undo the PW filter.
h0=GetNoiseWhitening(mi,nd,dis.ds);
h=h0;
h(h0<.25)=inf;
modelParticle=real(ifftn(fftn(filtParticle)./ifftshift(h)));

% Compute the residual image
resImg=ExtractImage(rscc.m0-rscc.mVes,[ix iy],nd)-modelParticle;

% Get the low-frequency values of the spectrum
spMask=Gaussian(nd,2,.1*nd); % pick up only the first few Fourier components
% (should be made resolution independent!)
me=(resImg(:)'*surroundMask(:))/sum(surroundMask(:)); % mean; do dc subtraction based on surround mask
mResImg=(resImg-me).*particleMask*dis.spectrumScale; % working residual image
fResImg=real(ifftn(fftn(mResImg).*ifftshift(h0))); % pre-whiten the residual
spc=fftshift(abs(fftn(fResImg)).^2).*spMask;
cropBoxSize=nd;
% spc=Crop(fftshift(abs(fftn(fResImg)).^2),cropBoxSize); % spectrum of cropped residual
cropBoxCtr=floor((cropBoxSize+1)/2);
spc(cropBoxCtr,cropBoxCtr)=0;
sval=max(spc(:))/(.2*partMaskRadius.^2); % The peak spectral estimate

% Get extrema values in the image
% We include 4 values which are the max(abs(gauss-filtered-image)) with
% filter frequencies of xfreqs (in A^-1), with and without residual
% subtraction.
% xvals  res subtr?  xfreqs
%    1      yes       .02
%    2      no        .02
%    3      yes       .01
%    4      no        .01
nxf=numel(xFreqs);
xvals=zeros(1,nxf);
rawImg=ExtractImage(rscc.m0-rscc.mVes,[ix iy],nd);
me=(rawImg(:)'*surroundMask(:))/sum(surroundMask(:)); % mean; do dc subtraction based on surround mask
mRawImg=(rawImg-me).*particleMask*dis.spectrumScale; % working residual image
ximgs=zeros(nd,nd,nxf,'single');
for j=1:2:nxf
    ximgs(:,:,j)=abs(GaussFilt(mRawImg,xFreqs(j)*pixA));
    ximgs(:,:,j+1)=abs(GaussFilt(mResImg,xFreqs(j)*pixA));
    xvals(j)=max(max(ximgs(:,:,j)))/(xFreqs(j)*10);
    xvals(j+1)=max(max(ximgs(:,:,j+1)))/(xFreqs(j)*10);
end;


if displayOn
    
    %     diagnostics on Fig. 4
    figure(4);
    subplot(221);
    imags(fResImg);
    subplot(222);
    imags(spc);
    subplot(223);
    plot(RadialPowerSpectrum(fResImg));
    
    
    figure(2);
    imul=2*dis.mulr;
    iadd=dis.addr;
    nd2=nd;
    nd2c=ceil((nd2+1)/2);
    br=dis.pars(20)/(mi.pixA*dis.ds*2);  % box radius
    boxX=br*[-1 -1 1 1 -1];
    boxY=br*[-1 1 1 -1 -1];
    
    %     Images for display
    mdImg=ExtractImage(rscc.m0,pos,nd2); % unsubtracted micrograph
    mdSub=mdImg-ExtractImage(rscc.mVes,pos,nd2); % subtracted micrograph
    partImg=Crop(modelParticle,nd2); % particle model
    resImg=mdSub-partImg;
    smsk=Crop(surroundMask,nd2);
    surroundImg=resImg.*smsk;
    me=sum(surroundImg(:))/sum(surroundMask(:));
    pmsk=Crop(particleMask,nd2);
    mResImg=(resImg-me).*pmsk;
    
    %     disp([s amp]);
    % Show the model particle
    fpartImg=GaussHP(GaussFilt(partImg,pixA/dis.filter(2)),pixA/dis.filter(1));
    fmdSub=GaussHP(GaussFilt(mdSub,pixA/dis.filter(2)),pixA/dis.filter(1));
    fmdImg=GaussHP(GaussFilt(mdImg,pixA/dis.filter(2)),pixA/dis.filter(1));
    
    mysubplot(361);
    imaga(fpartImg*imul+iadd); ShowBox;
    % Show the subtracted image, residual
    mysubplot(362); imaga(fmdSub*imul+iadd); ShowBox; axis off;
    %     title(amp);
    mysubplot(363); imaga((fmdSub-fpartImg)*imul+iadd); ShowBox; axis off;
    title(sval);
    %     Show the masked residual, original image, residual
    mysubplot(367); imaga((fmdSub-fpartImg).*pmsk*imul+iadd); ShowBox; axis off;
    mysubplot(368); imaga(fmdImg.*pmsk*imul+iadd); ShowBox; axis off;
    mysubplot(369); imaga((fmdImg-partImg).*pmsk*imul+iadd); ShowBox; axis off;
    
    % Close-up of the cc function
    mysubplot(3,6,13); imacs(ExtractImage(rscc.mxCCU,pos,nd/2)); axis off;
    mysubplot(3,6,14); imac(spc/(2*partMaskRadius)); axis off;
    rspc=Radial(spc)/(.2*partMaskRadius^2);
    mysubplot(3,6,15); bar(rspc); axis([1 cropBoxCtr-1 0 max(5,max(rspc))]);
    %     subplot(3,6,15); h=hist(spc(:)/128,0:2:80); bar(0:2:80,sqrt(h));
    mysubplot(122);
    imaga(GaussFilt(rscc.m0*imul+iadd,gfc));
    ShowBox(pos);
    axis off;
    title(mi.baseFilename,'interpreter','none');
    
    figure(5);
    for j=1:min(nxf,6)
        subplot(2,3,j);
        imags(ximgs(:,:,j));
        title(xvals(j));
    end;
    
    if displayOn>1
        pause
    end;
    figure(1);
end;

    function ShowBox(loc)
        if nargin<1
            loc=[nd2c nd2c];
        end;
        lw=2;
        hold on; plot(boxX+loc(1),boxY+loc(2),'y','linewidth',lw); hold off;
    end

% This function is copied from rsPickingPreprocessor4 so we can make the
% inverse filter.
    function T=GetNoiseWhitening(mi,np,ds)
        pwFiltPars=[.002 0; .01 .5];  % generic pw filter parameters
        fHP=.001; % highpass in A^-1
        T=meGetNoiseWhiteningFilter(mi,np,ds,1,fHP*mi.pixA*ds);
        if numel(T)<2
            %     disp('Using a generic filter');
            freqs=RadiusNorm(np)/pixA+1e-6;  % in inverse Å, prevent divide by zero
            T=ones(np,'single');
            for i=1:size(pwFiltPars,1)  % product of (gauss + const) terms.
                f0=pwFiltPars(i,1);
                a=pwFiltPars(i,2);
                h=exp(-(f0./freqs).^2);
                T=T.*(a+(1-a)*h);
            end;
        end;
        
    end



end
