% k2DamageCompNormalization
% Attempt to model the filtering--sum of all the filtered frames as in Grant & Grigorieff to correct the
% MotionCor2 output.
% Doesn't work...probably wrong.
% fs 18Jul17

ds=1;
segIndex=1;
startFrame=mi.frameSets(segIndex,1);
endFrame=mi.frameSets(segIndex,2);
nFrames=endFrame-startFrame+1;
n=mi.imageSize/ds;
f=RadiusNorm(n)/(mi.pixA*ds);  % 2d frequencies in A^-1
n0=min(1e5,eval(mi.damageModelCode));  % f is used in damageModelCode.
n2=2*n0;  % amplitude decay
k2=1./n2;  % decay rate
% Make the raw transfer functions as sqrt(SNR)
if numel(mi.frameDose)<2
    d=mi.frameDose;
else
    d=mean(mi.frameDose(startFrame:endFrame));
end;
d=.5;
% G=zeros(n,'single');
Qd=n2/d.*(1-exp(-d./n2))...
.*sqrt( (1-exp(-2*nFrames*d.*k2))./(1-exp(-2*d.*k2)) );  % sum of exposure transfer functions
Q=(1-exp(-nFrames*d*k2))./((1-exp(-d*k2)).*Qd);
% semilogy(sectr(Q));
% drawnow;

mcf=real(ifftn(fftn(mc).*ifftshift((1./Q))));
spf=RadialPowerSpectrum(mcf,0,4);
semilogy([spa spc spf]);
