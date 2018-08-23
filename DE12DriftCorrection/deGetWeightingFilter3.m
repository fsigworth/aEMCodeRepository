function [H1 rs]=deGetWeightingFilter3(wI1,wIMean,wISpec,hMask,wpixA,useCrystal)
% function [H rs]1=deGetWeightingFilter2(wI1,wIMean,wISpec,hMask,wpixA,useCrystal)
% From working-size images wI1 and corresponding independent means wIMean,
% create the weighting function for computing cross-correlations etc.hMask is a
% real-space mask applied to the images if desired, otherwise set to 1;
% however the correction with wISpec doesn't work with a mask.
% If the images contain a streptavidin crystal, the working pixel size and
% the useCrystal flag need to be set.
% We assume that the noise in the images is white. H1 is returned
% zero-orgin, and rs is a structure that gives the bandwidth-dependent
% varibles plotted in Shigematsu & Sigworth 2013, Fig. 4.

% Options
displayOn=0;
fch=.05;  % cross-spectrum gauss filter
fhp=.01;  % high-pass corner
mincsnr=2;
mincsnr2=2;

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
    trspect=zeros(nw);  % total real spectrum
    tvspect=zeros(nw);  % total imag power spectrum
    for i=1:nim
        im=double(wI1(:,:,i).*hMask);
        fMean=fftn(double(wIMean(:,:,i).*hMask));
        zspect=(fftn(im).*conj(fMean))/prod(nw)-wISpec(:,:,i);
        trspect=trspect+real(zspect);  % zero at origin
        tvspect=tvspect+2*imag(zspect).^2;
    end;
  
    ep=1e-6;
%     Blank the meridians where there is extra fixed noise.
    tvspect(1,:)=ep;
    tvspect(:,1)=ep;
    trspect(1,:)=ep;
    trspect(:,1)=ep;
%%
    csr1=Radial2(fftshift(trspect)/nim);  % radially averaged real amplitude
    csv1=Radial2(fftshift(tvspect))/nim;  % radially averaged pixel variance
rs.csr1=csr1;
rs.csv1=csv1;
    f=(1:numel(csr1))';
    w0=f.*csr1./csv1;   % initial weighting function W(f)
    w0(1)=0;
    filtW0=GaussFiltDCT(w0.*f,.02)./f;    % smoothed
    
    w0=max(w0,0);       % make it non-negative

%     Compute snr of the cc itself, ring by ring.
    rs.csnr=cumsum(w0.*csr1.*f).^2./cumsum(w0.^2.*csv1.*f);
    
    rs.csnr2=cumsum(w0.*csr1.*f.^3).^2./cumsum(w0.^2.*csv1.*f.^5);
    rs.csnr1=cumsum(w0.*csr1.*f.^3).^2./cumsum(w0.^2.*csv1.*f.^3);
    rs.w0=w0;

np=numel(w0);
% Limit the bandwidth, by two criteria.
% First, according to the 2nd derivative of the CC (sigma_2 in the paper)
% This is the criterion that seems to count.
p1=find(rs.csnr2>mincsnr2^2,1);  % find first point above threshold
p2=find(rs.csnr2>mincsnr2^2,1,'last'); % find last point above thresh.
if p2>p1+1  % nontrivial region above threshold
    w0(p2:np)=0;
end;

% Second, limit according to the CCF's SNR itself
p1=find(rs.csnr>mincsnr^2,1);  % find first point above threshold
p2=find(rs.csnr>mincsnr^2,1,'last'); % find last point above thresh.
if p2>p1+1
    w0(p2:np)=0;
end;

    
w2=Downsample(ToRect(w0),nw);
h0=(max(GaussFilt(w2,fch),0));

subplot(233);
% plot(sectr(h0));
mxa=max(sqrt([rs.csnr;rs.csnr2]));
s0=sectr(h0');
rs.s0=s0;
plot([sqrt(rs.csnr) sqrt(rs.csnr2) s0*mxa/max(s0)]);
title('Weighting filter');
drawnow;

    H1=fftshift(h0);
    H1(:,1)=0;  % blank the horizontal and vertical lines.
    H1(1,:)=0;
    norm=H1(:)'*H1(:)/numel(H1);
    H1=H1/sqrt(norm);  % rms mean value is 1.

end;
