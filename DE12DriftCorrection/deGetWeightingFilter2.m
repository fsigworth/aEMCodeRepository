function H1=deGetWeightingFilter2(wI1,wIMean,wISpec,hMask,wpixA,useCrystal)
% function H1=deGetWeightingFilter2(wI1,wIMean,wISpec,hMask,wpixA,useCrystal)
% From working-size images wI1 and corresponding independent means wIMean,
% create the weighting function for computing cross-correlations etc. as
% approximately amp(signal)/sd(noise), suppressing values <1.  wISpec is a
% stack of corresponding estimated stationary cross-spectra. hMask is a
% real-space mask applied to the images if desired, otherwise set to 1;
% however the correction with wISpec doesn't work with a mask
% If the images contain a streptavidin crystal, the working pixel size and
% the useCrystal flag need to be set.
% We assume that the noise in the images is white. H1 is returned
% zero-orgin.

% Options
maskLines=1;
displayOn=0;
fch=.05;  % cross-spectrum gauss filter
fhp=.01;  % high-pass corner

    wtConst=.05;

if nargin<4
    hMask=1;
end;
if nargin<5
    wpixA=0;
    useCrystal=0;
end;

[nx ny nim]=size(wI1);
nw=[nx ny];
hfMask=fuzzymask(nw,2,nw*.4,nw*.05);  % Final mask removing high frequencies

%% Compute the mean power spectrum
f=ifftshift(RadiusNorm(nw));  % frequencies, as singles
freqs=sectr(fftshift(f),2);  % 1D frequencies for radial spectrum. The second
% dimension is smaller than the first.
% Accumulate the 2d spectrum from all the images
sp2=zeros(nw);
for i=1:nim
    sp2=sp2+fftshift(abs(fftn(wI1(:,:,i))).^2);
end;
sp2=sp2/(nim*prod(nw));
sp1=Radial2(sp2);
nsp=numel(sp1);
spMean=mean(sp1(nsp/2:nsp));

if displayOn
    subplot(233);
    imacs(sp2.^.1);
    subplot(232)
    % sp=RadialPowerSpectrum(wI1(:,:,nim));
    loglog(freqs,[sp1/spMean sp1*0+1]);
    title('Spectrum');
end;

%% Make the weighting function
if useCrystal   % Just use the power spectrum, find the spots.
    spotSize=3;
    maxIndex=2;
    % wpixA=5.2;
    initialPars.a=52;
    initialPars.b=0;
    initialPars.alpha=pi/2;
    initialPars.theta=0;
    pars=deIndexLattice(sp2,wpixA,initialPars);
    spotPars=deGetSpotAmplitudes(sp2,wpixA,spotSize,pars,maxIndex);
    spotMask=deMakeSpotMask(nw,wpixA,spotSize,pars,spotPars);
    H1=ifftshift(sqrt(max(spotMask,0)));
    
else  % Create the weighting function from the cross-spectrum
    % First, accumulate the cross-spectrum beween images and means
    tcspect=zeros(nw);
    for i=1:nim
        im=double(wI1(:,:,i).*hMask);
        fMean=fftn(double(wIMean(:,:,i).*hMask));
        cspect=real(fftn(im).*conj(fMean))/prod(nw);
        tcspect=tcspect+cspect;  % zero at origin
    end;
    tcspect=fftshift(real(tcspect-sum(wISpec,3)));  % Remove the mean excess noise.
    % now is zero at center
    if maskLines  % Get rid of vertical and horizontal lines
        tcspect(nx/2+1,:)=0;
        tcspect(:,ny/2+1)=0;
    end;
    csp1=Radial2(tcspect);   % get the 1d function
    ncs=numel(csp1);
    csp1=csp1-mean(round(ncs*.9):ncs);
    
    snr1=GaussFilt(csp1./sp1,.05);
    snr1=max(snr1,1e-3);
    subplot(233);
    semilogy(snr1);
    title('Radial SNR');
    drawnow;
    
    
    msk0=ToRect(exp(-wtConst./sqrt(snr1)));  %% This is the weighting function
    msk=Downsample(msk0,nw).*(fftshift(f));
    %
    %     ftry=0.05;
    %     minSNR=.02;
    %     maxf=find(diff(GaussFiltDCT(snr1,ftry)>minSNR)~=0,1,'first');
    %     nfound=numel(maxf);
    % %     while nfound>1 && ftry > .01
    % %         ftry=ftry/sqrt(2);
    % %         maxf=find(diff(GaussFilt(snr1,ftry)>minSNR)~=0);
    % %         nfound=numel(maxf);
    % %     end;
    % %     maxf=maxf(1);
    %     if numel(maxf)<1  % no points below 0.1!
    %         maxf=ny/2
    %     elseif maxf<ny/4  % at least allow frequencies up to .1/dsw
    %         maxf=ny/4
    %     end;
    %     maxf
    % maxf=nw*0.4;
    %     msk=fuzzymask(nw,2,maxf,maxf/2);
    
    % High-pass
    minf=fhp*nw;
    msk=msk.*(1-fuzzymask(nw,2,minf,minf/10));
    msk=msk.*hfMask;
    %     msk=msk.*fftshift(f);   % include some f dependence.
    %     tcspect=tcspect.*msk;   % zero is at origin
    tcspect=msk;
    %% compute the masking function
    % Smooth the 2D spectrum
    h0=(max(GaussFilt(tcspect,fch),0))/(spMean);
    H1=fftshift(h0);
    H1(:,1)=0;  % blank the horizontal and vertical lines.
    H1(1,:)=0;
    norm=H1(:)'*H1(:)/numel(H1);
    H1=H1/sqrt(norm);  % rms mean value is 1.
    
end;
