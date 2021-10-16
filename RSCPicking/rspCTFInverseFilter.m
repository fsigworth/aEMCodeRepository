function m1=rspCTFInverseFilter(m,mi,compFraction,effWeights)
% function m=rspCTFInverseFilter(m,mi,comp)
% Given the downsampled image m, boost the low frequencies below the 1st
% peak of the first exposure's CTF.  If compFraction=1, the boost is the inverse of
% the ctf, up to a value of 10.
% If we're padding the image m, it would be good to have it correspond to
% mi.padImageSize.
% Includes code to make a fake double-exposure image if weights are not all 1.
%
% mi fields that are used
%   weights
%   pixA
%   padImageSize (or imageSize if padImageSize not given)
%   ctf
if nargin<5 || sum(effWeights)==0
    effWeights=ones(size(mi.weights));
end;
doDoubleExposure = ~all(effWeights(:)==mi.weights(:)); % emulate double exposure
% if some of mi.weights are not =1.

if compFraction==0 && ~doDoubleExposure % nothing to do
    m1=m;
    return
end;

k=.1;

n=size(m);
if ~isfield(mi,'padImageSize');
    origSize=mi.padImageSize;
else
    origSize=mi.imageSize;
end;    
    ds=origSize(1)/n(1);  % actual downsampling of m
pixA=mi.pixA*ds;

freqs=RadiusNorm(n)/pixA;  % frequencies of padded image
[ctf,chi]=ContrastTransfer(freqs,mi.ctf(1));  % Just get chi
peakMask=abs(chi)<pi/4;  % mask up to just beyond the first peak
mxAmp=max2d(peakMask.*abs(ctf));
H=ones(n,'single');
if compFraction>0
    H(peakMask)=mxAmp*(1-compFraction+compFraction./(abs(ctf(peakMask))+k));
end;
% Hp=GaussHPKernel(n,fHP);

%%%%%%
if doDoubleExposure
    c0=meGetEffectiveCTF(mi,n,ds);
    mi.weights=effWeights;
    c1=meGetEffectiveCTF(mi,n,ds);
    H=H.*c1./c0;
%    H=H.*(1+weightComp*(c1./c0-1));
% figure(2); clf; plot([sectr(c1) sectr(c0) sectr(c1./c0)]); figure(1);

end;
%%%%%%%%%

m1=real(ifftn(fftn(m).*ifftshift(H)));
