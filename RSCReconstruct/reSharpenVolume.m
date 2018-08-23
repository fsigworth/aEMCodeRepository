% reSharpenVolume.m
% Sharpen the Kv volume tVol.  We need pixA and ri.volMask.

% we read the Kv1.2 map, and compare its spectrum within the volume mask
% with the tVol map.

%%
si.pixA=2.4;
pixA=si.pixA;
% tVol=1.67*Downsample(Crop(vVol,64),128)+Crop(pVol,128);

% tVol=pVol;

n=size(tVol,1);
freqs=(1:n/2)'/(n*pixA);
msk3=fuzzymask(n,3,[.18 .18 .12]*n,0.04*n,[0 0 -0.12*n]+(n/2+1));
figure(1);
ShowSections(tVol.*msk3,[],45);
tSpect=RadialPowerSpectrum(tVol.*msk3);
mysubplot(3,3,9);
semilogy(freqs,tSpect);
drawnow
%%
kMag=.92;
q=load('/Users/fred/Structures/kv1.2/KvMap.mat');
kvMap=DownsampleGeneral(q.sim,n,kMag*q.pixA/si.pixA);
rKvMap=ERotate3(kvMap,[0 pi/2 pi/2]);
figure(2);
ShowSections(rKvMap.*msk3,[],45);
rSpect=RadialPowerSpectrum(rKvMap.*msk3);
mysubplot(3,3,9);
semilogy(freqs,rSpect);
drawnow;
%%
df=1/(n*pixA);
figure(3);
f0=floor(.03/df);
tNorm=mean(tSpect(1:ceil(f0)));
rNorm=mean(rSpect(1:ceil(f0)));

rSpect=rSpect*tNorm/rNorm;
mysubplot(221);
semilogy(freqs,[tSpect rSpect]);
fc=.09;
fw=.02;
fe=.015;
fl=.06;
filt2=rSpect./tSpect.*(0.5+0.5*erf((fc-freqs)/fw));
filt=exp((max(0,freqs-fl))/fe).*(0.5+0.5*erf((fc-freqs)/fw));
filt(1:f0)=1;
mysubplot(222);
plot(freqs,[filt filt2]);

%%
R3=Radius3(n);
F3=filt2(max(1,min(n/2,round(R3))));

fVol=real(ifftn(fftn(tVol).*ifftshift(F3)));
xMsk=Crop(ri.volMask,n);
% xMsk=msk3;
% xMsk=1;
figure(4);

ShowSections(fVol.*xMsk,[],45)




doWrite=1;

if doWrite
WriteMRC((fVol.*xMsk),2.4,'R112k6v0MbnSharp02Mask.mrc');
% WriteMRC((tVol.*xMsk),2.4,'R112k6v0MbnMask.mrc');
end;

return
%%
WriteMRC(fVol,si.pixA,'i29av02sharp.mrc');
