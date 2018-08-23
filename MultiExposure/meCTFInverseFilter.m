function [mfilt,H]=meCTFInverseFilter(m,mi,lfAmp,fDeTrend,fHP)
% function mfilt=meCTFInverseFilter(m,mi,lfAmp,fDeTrend,fHP)
% Flatten the overall transfer function below the first zero of the
% effective CTF. lfAmp is the final boost of LF components--nominally 1,
% but can be reduced to reveal HP components better.
% fDeTrend is the frequency (A^-1) of the de-trending filter that avoids
% edge artifacts; default is .0005.  To preserve DC information (e.g.
% density of carbon in image) set this to zero.
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

padFactor=1.5;
n=size(m);
n1=n*padFactor;  % scale up to avoid edge artifacts
ds=mi.imageSize(1)/n(1);  % actual downsampling of m
pixA=mi.pixA*ds;

if fDeTrend>0
    mdc=GaussFiltDCT(m,fDeTrend*pixA);
    m1=m-mdc;
else
    m1=m;
end;

freqs=RadiusNorm(n1)/pixA;  % frequencies of padded image
[ctf,chi]=ContrastTransfer(freqs,mi.ctf(1));  % Just get chi to find first zero.
effCTF=meGetEffectiveCTF(mi,n1,ds);
peakMask=abs(chi)<.5;  % frequencies below the first maximum.
peakPtr=floor(n1(1)/2)+find(-sectr(chi)>.5,1);
if numel(peakPtr)<1
    peakPtr=1;
end;
zeroMask=abs(chi)<1;   % below the first zero
highPass=n1(1)*pixA*fHP;
H=ones(n1);
% Set a constant value between the peak and the first zero
H(zeroMask)=lfAmp/effCTF(peakPtr,ceil((n1(2)+1)/2));
% Reciprocal of effCTF below the peak frequency
H(peakMask)=lfAmp./effCTF(peakMask);
H=H.*(1-fuzzymask(n1,2,highPass,highPass/4));
mfilt=Crop(real(ifftn(fftn(Crop(m1,n1)).*ifftshift(H))),n);
%     imacs(mfilt);

