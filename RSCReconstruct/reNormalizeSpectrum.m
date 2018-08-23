% reNormalizeSpectrum
% Sharpen 3d reconstructions, do 3D alignments and create 3D masks
% For Katrine's JSB paper

fc=0.9;
zSlice=54;


% cd('~/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/Recon96jl/mrc')
[m,d]=ReadMRC('i29cv01.mrc');
n=size(m,1);
ctr=ceil((n+1)/2)
fc1=fc*d.pixA;
% m=SharpFilt(m,fc1,fc1/10);

figure(1);
dp1=RadialPowerSpectrum(m,0);
fs=(1:n/2)'/(96*d.pixA);
figure(1);
ShowSections(m,[],45);
subplot(3,3,9);
semilogy(fs,dp1)
drawnow;

s=load('~/aEMCodeRepository/AMPAR/KvMap.mat');
ns=size(s.map,1);
fc1=fc*d.pixA;
refMap=DownsampleGeneral(s.map,n,s.pixA/d.pixA);
% refMap=SharpFilt(refMap,fc1,fc1/10);

sp1=RadialPowerSpectrum(refMap,0);

refMapAli=reAlignVolumes(m,refMap);
figure(2);
ShowSections(SharpFilt(refMapAli,fc1),[ctr ctr zSlice],45);
subplot(3,3,9);
semilogy(fs,sp1);

%%
msk=GaussFilt(GaussFilt(refMapAli>.1,.1)>.2,.1);
mmsk=m.*msk;
% mmsk=SharpFilt(m.*msk,fc1);

dm1=RadialPowerSpectrum(mmsk,0);

figure(3);
ShowSections(mmsk,[],45);
subplot(3,3,9);
semilogy(fs,dm1);

sscl=sp1(1)/dm1(1);

% p=ceil(fc1*n);
% dm1(p:end)=dm1(p);
% sp1(p:end)=sp1(p);

figure(4);
semilogy(fs,[dm1*sscl sp1]);

%%
% make the whitening filter
w1=sqrt(sp1./(dm1*sscl));
% w3=zeros([n n n],'single');
R=round(Radius3(n));
R(R<1)=1;
R(R>n/2)=n/2;
w3=w1(R);
w3=GaussFilt(w3,.05);
ShowSections(w3);
subplot(3,3,9);
plot(fs,Radial3(w3));

%% Apply it
mf=real(ifftn(fftn(m).*ifftshift(w3)));
mfmsk=mf.*msk;
mf2msk=SharpFilt(mfmsk,fc1);
mf2=SharpFilt(mf,fc1);
figure(3);
ShowSections(mf2msk,[ctr ctr zSlice],45);
% ShowSections(mf2,[],45);
dmf=RadialPowerSpectrum(mf2msk,0);
subplot(3,3,9);
semilogy(fs,[dmf*sscl sp1]);


%% Apply it to membrane-restored model
msk2=GaussFilt(GaussFilt(refMapAli>.1,.05)>.1,.1);

WriteMRC(mf2.*msk2,d.pixA,'i29cv01wmsk.mrc');
WriteMRC(mf2,d.pixA,'i29cv01w.mrc');
mv=ReadMRC('i30cv01.mrc');
mvf=real(ifftn(fftn(mv).*ifftshift(w3)));
mvf2=SharpFilt(mvf,fc1);
ShowSections(mvf2.*msk2,[ctr ctr zSlice],45);

WriteMRC(mvf2.*msk2,d.pixA,'i30cv01wmsk.mrc');
WriteMRC(mvf2,d.pixA,'i30cv01w.mrc');
WriteMRC(msk2,d.pixA,'i30cv01mask.mrc');
