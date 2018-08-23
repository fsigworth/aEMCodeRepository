% reCrystalFSC.m


xMapName='~/Structures/kv1.2/KvMapMbnSub.mrc';
eMapPath='/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Reconstructions/';
eMapExpt='Recon96imi/mrc/';
eMapExpt='Recon96ir/mrc/';
eMapName='i29cv01.mrc';
eMapName='i30cv01.mrc';



[xMap0,s]=ReadMRC(xMapName);
xPixA=s.pixA;

[eMap0,s]=ReadMRC([eMapPath eMapExpt eMapName]);
pixA=s.pixA;
n=size(eMap0,1);

mag=xPixA/pixA;
mags=mag*[.95 .975 1 1.025];
mags=mag*.98;
fscs=zeros(n/2,numel(mags)+1);

for ind=1:numel(mags)
    
xMap=DownsampleGeneral(shiftdim(xMap0,2),n,mags(ind));

msk=fuzzymask(n,3,.35*n,.2*n);

[eMap,tz,gamma,mirror]=reAlignVolumes(xMap,msk.*eMap0,[0 0 1 1]);
tz
gamma
mirror

msk2=fuzzymask(n,3,.2*n,.1*n, (n/2+1)+[0 0 -.1*n]);
msk2=1;

figure(1);
ShowSections(eMap.*msk2,[],45);
figure(2);
ShowSections(xMap,[],45);


fsc=FSCorr2(eMap.*msk2,xMap);
fsc(1)=1;
fsc=(fsc+circshift(fsc,[1 0])+circshift(fsc,[-1 0]))/3;
fsc(1)=1;
fscs(:,ind)=fsc;
figure(1);
subplot(339);
freqs=(.5:n/2-.5)/(n*pixA);
plot(freqs,fscs);
axis([0 .15 -.05 1]);
drawnow;
end;