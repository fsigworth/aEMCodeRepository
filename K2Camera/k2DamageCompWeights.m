function [H,effCTF]=k2DamageCompWeights(mi,segIndex,ds)
% function [H,effCTF]=k2DamageCompWeights(mi,segIndex,ds)
% segIndex is the index into rows of mi.frameSets.  ds is the optional
% downsampling factor, if you are working on downsampled frames.
% H is the n x n x nFrames stack of transfer functions, zero in center.
% Compute the transfer functions for each frame of a movie segment.
% Form the accumulated image by scaling the ith frame by H(:,:,i), e.g.
%   fo=real(ifft2(fft2(fi).*ifftshift2(H)));, then summing.
% or, probably better to use doubles to avoid roundoff errors:
%   fo=single(real(ifft2(fft2(double(fi))).*double(ifftshift2(H))));
% Changed 28 Oct 14 to allow mi.frameDose to be a vector.
% Changed 30 Oct 16 to use 1/(2*n0) instead 1/n0 decay of amplitudes!

% % test code
% mi=meCreateMicrographInfoStruct14;
% mi.imageSize=[3840 3840];
% mi.pixA=1.247;
% mi.frameSets=[2 21; 22 40];
% mi.frameDose=1.5;
% for segIndex=1:2
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
n0=min(1e5,eval(mi.damageModelCode));  % f is used in damageModelCode.
n2=2*n0;  % amplitude decay
k2=1./n2;  % decay rate
% Make the raw transfer functions as sqrt(SNR)
if numel(mi.frameDose)<2
    d=mi.frameDose*ones(nFrames,1);
else
    d=mi.frameDose(startFrame:startFrame+nFrames-1);
end;
accumDose=[0 cumsum(d(:))'];
H=zeros([n nFrames],'single');
% G=zeros(n,'single');
for i=1:nFrames
    Q=n2/d(i).*(1-exp(-d(i)*k2));  % average decay during the exposure
    H(:,:,i)=Q.*exp(-accumDose(i)*k2);
%     H(:,:,i)=exp(-accumDose(i)*k2).*Q;
end;
% normalize to give white shot noise
noiseGain=sqrt(sum(H.^2,3)/nFrames);  % scaling of white noise spectrum
effCTF=noiseGain;
for i=1:nFrames
    H(:,:,i)=H(:,:,i)./noiseGain;
end;

%%
% % More test code
% uncorrCTF=n0/(nFrames*d).*exp(-d*(startFrame-1)*k0).*(1-exp(-(d*nFrames*k0)));
% subplot(3,2,segIndex)
% plot(sectr(f),[sectr(effCTF) sectr(uncorrCTF)]*sqrt(nFrames));
% axis([0 .4 0 inf]);
% legend('With damage comp','Without');
% xlabel('Spatial frequency, A^{-1}');
% ylabel('Effective CTF');
% title(['Frames ' num2str(mi.frameSets(segIndex,:))]);
% subplot(3,2,segIndex+2);
% plot(sectr(f),sectr(H));
% ylabel('Weights');
% axis([0 .4 0 inf]);
% 
% subplot(3,2,segIndex+4);
% plot(sectr(f),sectr(n0));
% ylabel('Critial dose (amplitude)');
% axis([0 .4 0 inf]);
% 
% end;