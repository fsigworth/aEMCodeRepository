function [mfilt,H]=meCTFInverseFilter(m,mi,lfAmp,fDeTrend,fHP)
% function mfilt=meCTFInverseFilter(m,mi,lfAmp,fDeTrend,fHP)
% Flatten the overall transfer function below the first maximum of the
% effective CTF. lfAmp is the final boost of LF components. lfAmp=1 means
% a complete inverse filter from zero to the first peak, then unity at all
% higher frequencies.
% m is assumed to be derived from a padded image, if mi.padImageSize is
% present.
% fDeTrend is the frequency (A^-1) of the de-trending filter that avoids
% edge artifacts; default is .0005.  To preserve DC information (e.g.
% density of carbon in image) set this to zero. H does not include the
% detrend filter.
% fHP is the high-pass corner (A^-1) of
% the inverse filter; default is also .0005.
% The (zero centered) filter transfer function is returned as H.

if nargin<3
    lfAmp=1;
end;
if nargin<4
    fDeTrend=.0005;
end;
if nargin<5
    fHP=.0005;
end;

useSimpleCTF=1; % Don't compute the full effective CTF.
padFactor=1.5;
n=size(m);
n1=n*padFactor;  % we will pad in Fourier domain to avoid edge artifacts
imSize=mi.imageSize;
if isfield(mi,'padImageSize')
    imSize=mi.padImageSize;
end;
ds=imSize(1)/n(1);  % actual downsampling of m

pixA=mi.pixA*ds; % pixel size of m

if fDeTrend>0
    mdc=GaussFiltDCT(m,fDeTrend*pixA);
    m1=m-mdc;
else
    m1=m;
end;

freqs=RadiusNorm(n1)/pixA;  % frequencies of padded image

[effCTF,chi]=ContrastTransfer(freqs,mi.ctf(1));  % Just get chi to find first peak.
if ~useSimpleCTF
    effCTF=meGetEffectiveCTF(mi,n1,ds);
end;

peakMask=-chi<.5;  % frequencies below the first maximum.
peakPtr=floor(n1(1)/2)+find(-sectr(chi)>.5,1); % a peak point on the x axis
if numel(peakPtr)<1
    peakPtr=1;
end;
peakVal=effCTF(peakPtr,ceil((n1(2)+1)/2)); % a value close to 1
H=ones(n1);
% Reciprocal of effCTF below the peak frequency
excessInverse=peakVal./effCTF(peakMask)-1;
H(peakMask)=lfAmp*excessInverse+1; % linear scaling the boost beyond 1

highPass=n1(1)*pixA*fHP;
H=H.*(1-fuzzymask(n1,2,highPass,highPass/4)); % insert the high pass
mfilt=Crop(real(ifftn(fftn(Crop(m1,n1)).*ifftshift(H))),n);
%     imacs(mfilt);

