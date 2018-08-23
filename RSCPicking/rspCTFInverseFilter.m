function m1=rspCTFInverseFilter(m,mi,comp,effWeights)
% function m=rspCTFInverseFilter(m,mi,comp)
% Given the downsampled image m, boost the low frequencies below the 1st
% peak of the first exposure's CTF.  If comp=1, the boost is the inverse of
% the ctf, up to a value of 10.
%%%%%% include fake double-exposure code.
if nargin<5
    effWeights=ones(size(mi.weights));
end;
doDoubleExposure = ~all(effWeights(:)==mi.weights(:)); % emulate double exposure
if comp==0 && ~doDoubleExposure % nothing to do
    m1=m;
    return
end;

k=.1;

n=size(m);
ds=mi.imageSize(1)/n(1);  % actual downsampling of m
pixA=mi.pixA*ds;

freqs=RadiusNorm(n)/pixA;  % frequencies of padded image
[ctf,chi]=ContrastTransfer(freqs,mi.ctf(1));  % Just get chi
peakMask=abs(chi)<pi/4;  % mask up to just beyond the first peak
mxAmp=max2d(peakMask.*abs(ctf));
H=ones(n,'single');
if comp>0
    H(peakMask)=mxAmp*(1-comp+comp./abs(ctf(peakMask))+k);
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
