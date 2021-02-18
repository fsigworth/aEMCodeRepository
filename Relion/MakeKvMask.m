% MakeKvMask.m
outputName='KvMask128.mrc';
n=128;
n0=128; %??
pixA=2.26
us=1.1;
% pa='/Users/fred/Structures/kv1.2/';
pa='~/data/';
% [m0,s]=ReadMRC([pa 'Kv1.2Tetramer64.mrc']);
[m0,s]=ReadMRC([pa 'KvRef1.05.mrc']);
% m0=shiftdim(m0,2);
% m0=circshift(m0,[0 0 3]);
ms=DownsampleGeneral(m0,n0,s.pixA/pixA*us);
%% Not very good mask here...
% msk1a=GaussFilt((ms),.1)>.1;
msk1b=GaussFilt((ms),.05)>.01;
% tmb=round(0.55*n);
% msk1a(:,:,tmb:end)=msk1b(:,:,tmb:end);
msk=GaussFilt(msk1b,.1);
figure(2);
ShowSections(msk.*ms,[],45);
figure(3);
ShowSections(ms,[],45);
figure(4);
ShowSections(msk,[],45);
% WriteMRC(msk,s.pixA,outputName);
