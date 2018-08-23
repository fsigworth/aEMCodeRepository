% MakeKvMask.m
outputName='KvMask.mrc';
n=96;
pixA=3.066;
us=1.1;
[m0,s]=ReadMRC('Kv1.2Tetramer64.mrc');
m0=shiftdim(m0,2);
m0=circshift(m0,[0 0 3]);
ms=DownsampleGeneral(m0,n0,s.pixA/pixA*us);
%%
msk1a=GaussFilt((ms),.1)>.1;
msk1b=GaussFilt((ms),.1)>.01;
tmb=round(0.55*n);
msk1a(:,:,tmb:end)=msk1b(:,:,tmb:end);
msk=GaussFilt(msk1a,.1);
figure(2);
ShowSections(msk.*ms,[49 49 57],45);
figure(3);
ShowSections(ms,[49 49 57],45);
figure(4);
ShowSections(msk,[49 49 57],45);
WriteMRC(msk,s.pixA,outputName);
