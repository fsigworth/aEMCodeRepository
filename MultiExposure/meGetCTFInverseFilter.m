function h=meGetCTFInverseFilter(mi,n,ds,lfGain)
% function h=meGetCTFInverseFilter(mi,n,ds,hfGain)
% Compute a filter transfer function that will flatten the low-frequency
% contrast-transfer function of a merged image. The transfer function h
% below the first peak is the inverse filter 
% lfGain * effCTF(1stPeak)/effCTF(f); between the first peak and first zero
% it is lfGain; at high frequencies it is unity.
k=0;

if numel(n)~=2
    n=[n(1) n(1)];
end;
if nargin<3
    ds=mi.imageSize(1)/n(1);
end;
if nargin<4
    lfGain=1;
end;
if ~isfield(mi,'weights');
    mi.weights=ones(size(mi.doses));
end;

effPixA=mi.pixA*ds;
freqs=RadiusNorm(n)/effPixA;  % frequencies for evaluating the CTF
% Find the first zero of the first effective exposure
iexp1=find(mi.weights.*mi.doses,1,'first');
[ctf,chi]=ContrastTransfer(freqs,mi.ctf(iexp1));

effCTF=meGetEffectiveCTF(mi,n,ds);

fBelowPeak=chi>-.5;  % frequency mask beyond first max.
fBelowZero=chi>-1;  % frequencies beyond first zero.
[q, ipk]=min((chi(:)+.5).^2);
apk=effCTF(ipk);
h=ones(n);
h(fBelowZero)=lfGain;
h(fBelowPeak)=lfGain*(apk+k)./(effCTF(fBelowPeak)+k);
