
cd('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon112j2/mrc')

[v1,s]=ReadMRC('i68av01.mrc');
v2=ReadMRC('i68bv01.mrc');
ShowSections(v1)
n=size(v1,1);
ctr=ceil((n+1)/2);
volMask=min(1,fuzzymask(n,3,[.32 .32 .15]*n,.1*n,[ctr ctr ctr+.18*n])...
                    +fuzzymask(n,3,0.22*n,0.1*n,[ctr ctr ctr-.06*n]));
ShowSections(v1.*volMask);
v1m=v1.*volMask;
v2m=v2.*volMask;


[v2ali,tz,gamma,mirror]=reAlignVolumes(v1m,v2);
disp([tz gamma mirror]);

%%
v0=(v1m+v2ali.*volMask)/2;
vu=(v1+v2ali)/2;

% WriteMRC(v0,s.pixA,'i68dv01.mrc')
%%
freqs=(0:(n/2)-1)/(n*s.pixA);
fsc=FSCorr2(v2ali.*volMask,v1m);
figure;
plot(freqs,fsc);

% normalize the power spectrum
%%
load('/Users/fred/aEMCodeRepository/AMPAR/KvMap.mat');
q=1.2;
refMap=DownsampleGeneral(map,n,q*pixA/s.pixA);
ShowSections(refMap);
%%
s=RadialPowerSpectrum(v0);
s0=RadialPowerSpectrum(refMap);
s0=s0*s(2)/s0(2);
semilogy(freqs,[s s0]);
%%
f0=.09;
fw=.02;
fmsk=.5-.5*erf((freqs'-f0)/fw);
plot(freqs,fmsk);

w=sqrt(s0.*fmsk./s);
plot(freqs,w);

R=Radius3(n);
W=interp1(w,R(:));
W(isnan(W))=0;
W=reshape(W,n,n,n);

v0f=real(ifftn(fftn(v0).*ifftshift(W)));
vuf=real(ifftn(fftn(vu).*ifftshift(W)));
figure(1);
ShowSections(vuf);

%%
% WriteMRC(v0f,2.5,'i68w09.mrc');
%%
figure(5);
ImagicDisplay3(vuf);

