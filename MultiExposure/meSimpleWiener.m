function [mcf h]=meSimpleWiener( mc, effctz, pixAd, fc, displayOn )
% function out=meSimpleWiener( mc, effctz, pixAd, fc, displayOn )
% Correct the effective ctf down to low frequencies for a merged image mc.
% effctz is the effective CTF as returned by meComputeMergeCoeffs and pixAd
% is the pixel size of mc.
% fc is the frequency of transition from wiener (at low frequencies) to no
% filter (at high frequencies).  An example might be fc=.015 A^-1.
% fc=.015;
if nargin<5
    displayOn=0;
end;
fcd=fc/2;  % spread of cutoff
[nd ny]=size(mc);
df=1/(nd*pixAd);  % frequency increment
swf=fuzzymask(nd,2,fc/df,fcd/df);  % switch function
gsw=-diff(sectr(swf));
gsw(nd/2)=0;

MeanGain=sum(gsw.*sectr(effctz))/sum(gsw);
epsi=.001;
h=MeanGain*effctz./(effctz.^2+epsi).*swf+(1-swf);

mcf=real(ifftn(fftn(mc).*fftshift(h)));


if displayOn  % Plot the various transfer functions
    f=(0:nd/2-1)'*df;
    semilogx(f,[sectr(effctz) sectr(h), sectr(h.*effctz)]);
    legend('Merged eff CTF', 'Filter', 'Result');
    xlabel('Frequency, A^{-1}');
    ylabel('Transfer function');
end;
%
% figure(2); SetGrayscale;
%
% img=imscale(GaussFilt(-mc,.2));
% imac((img-100)*2);