% MakeKvMask.m
outputName='EMWork/KvMask192pixA2.mrc';
looseOutputName='EMWork/KvMask192pix2Loose.mrc';
n=192;
n0=192; %??
pixA=2;
us=1;
pa='/Users/fred/Structures/kv1.2/';
[m0,s]=ReadMRC([pa 'Kv1.2Tetramer64.mrc']);
m0=shiftdim(m0,2);
m0=circshift(m0,[0 0 3]);
disp(['Map pixel size is ' num2str(s.pixA)]);
ms=DownsampleGeneral(m0,n0,s.pixA/pixA*us);
%%
msk1a=GaussFilt((ms),.1)>.1;
msk1b=GaussFilt((ms),.1)>.01;
tmb=round(0.55*n);
msk1a(:,:,tmb:end)=msk1b(:,:,tmb:end);
msk=GaussFilt(msk1a,.05);
figure(2);
% ShowSections(msk.*ms,[49 49 57],45);
% figure(3);
% ShowSections(ms,[49 49 57],45);
% figure(4);
% ShowSections(msk,[49 49 57],45);
msk=min(1,max(0,msk)); % clip for Relion
WriteMRC(msk,pixA,outputName);
%%
blurMsk=GaussFilt(msk,.1)>.01;
blurMsk=min(1,max(0,blurMsk)); % clip for Relion
ShowSections(blurMsk,[],45);
WriteMRC(blurMsk,pixA,looseOutputName)