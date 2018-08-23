function effCTF=k2DamageCompCTF(mi,segIndex,ds)
% function effCTF=k2DamageCompCTF(mi,segIndex,ds)
% segIndex is the index into rows of mi.frameSets.  ds is the optional
% downsampling factor, if you are working on downsampled frames.
% effCTF is the zero-center CTF factor from damage after damage
% compensation.

% % test code
% mi=meCreateMicrographInfoStruct14;
% mi.imageSize=[3840 3840];
% mi.pixA=1.247;
% mi.frameSets=[1 20; 22 40];
% mi.frameDose=1.5;
% segIndex=1;
% %

if nargin<3
    ds=1;
end;
if mi.version<14
    error(['mi structure version should be >=14, is ' num2str(mi.version)]);
end;

startFrame=mi.frameSets(segIndex,1);
nFrames=mi.frameSets(segIndex,2)-startFrame+1;
n=mi.imageSize/ds;
f=RadiusNorm(n)/(mi.pixA*ds);  % 2d frequencies in A^-1
n0=eval(mi.damageModelCode);
k0=1./n0;  % decay rate
% Make the raw transfer functions as sqrt(SNR)
d=mi.frameDose;
G=zeros(n,'single');
Q=1-exp(-d*k0);
for i=1:nFrames
    dStart=d*(startFrame+i-2);
    h=(n0/d).*exp(-dStart*k0).*Q;
    G=G+h.^2;
end;
% normalize to give white shot noise
effCTF=sqrt(G/nFrames);



% % More test code
% uncorrCTF=n0/(nFrames*d).*exp(-d*(startFrame-1)*k0).*(1-exp(-(d*nFrames*k0)));
% subplot(2,1,1)
% plot(sectr(f),[sectr(effCTF) sectr(uncorrCTF)]);
% axis([0 inf 0 1]);
% legend('With damage comp','Without');
% xlabel('Spatial frequency, A^{-1}');
% ylabel('Effective CTF');
% 
% subplot(2,1,2);
% plot(sectr(f),sectr(n0));
