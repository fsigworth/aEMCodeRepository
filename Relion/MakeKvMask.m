% MakeKvMask.m
% outputName='EMWork/KvMask192pixA2.mrc';
% looseOutputName='EMWork/KvMask192pix2Loose.mrc';
% n=192;
% n0=192; %??
% pixA=2;
% us=1;
% pa='/Users/fred/Structures/kv1.2/';
% [m0,s]=ReadMRC([pa 'Kv1.2Tetramer64.mrc']);
% m0=shiftdim(m0,2);
% m0=circshift(m0,[0 0 3]);
% disp(['Map pixel size is ' num2str(s.pixA)]);
% ms=DownsampleGeneral(m0,n0,s.pixA/pixA*us);
% %%
% msk1a=GaussFilt((ms),.1)>.1;
% msk1b=GaussFilt((ms),.1)>.01;
% tmb=round(0.55*n);
% msk1a(:,:,tmb:end)=msk1b(:,:,tmb:end);
% msk=GaussFilt(msk1a,.05);
% figure(2);
% % ShowSections(msk.*ms,[49 49 57],45);
% % figure(3);
% % ShowSections(ms,[49 49 57],45);
% % figure(4);
% % ShowSections(msk,[49 49 57],45);
% msk=min(1,max(0,msk)); % clip for Relion
% WriteMRC(msk,pixA,outputName);
% %%
% blurMsk=GaussFilt(msk,.1)>.01;
% blurMsk=min(1,max(0,blurMsk)); % clip for Relion
% ShowSections(blurMsk,[],45);
% WriteMRC(blurMsk,pixA,looseOutputName)

tightMaskName='';
looseMaskName='';
ellipseMaskName='';
discMaskName='KvMaskCirc.mrc';
outPath='RSC9/'

%%%% farnam version

inputName='KvRef128.mrc';
% pa='/Users/fred/Structures/kv1.2/'; % Katz
inputPath='~/mini_drobo/' % siggpu

n=128
pixA=2.26
us=1.1;
% [m0,s]=ReadMRC([pa 'Kv1.2Tetramer64.mrc']);
[m0,s]=ReadMRC([inputPath 'KvRef.mrc']);
% m0=shiftdim(m0,2);
% m0=circshift(m0,[0 0 3]);
ms=DownsampleGeneral(m0,n0,s.pixA/pixA*us);
%% Tight mask
msk1b=GaussFilt((ms),.05)>.01;
msk=GaussFilt(msk1b,.1);
figure(2);

% ShowSections(msk.*ms,[49 49 57],45);
msk=min(1,max(0,msk)); % clip for Relion
if numel(tightMaskName)>0
    WriteMRC(msk,pixA,tightMaskName);
    disp([tightMaskName ' written.']);
end;
%% Looser mask
blurMsk=GaussFilt(GaussFilt(msk,.1)>.01,.1);
blurMsk=min(1,max(0,blurMsk)); % clip for Relion
ShowSections(blurMsk,[],45);
% if numel(looseMaskName)>0
%     WriteMRC(blurMsk,pixA,looseMaskName)
% end;

%% With membrane ellipsoid
sphere=fuzzymask(n,3,n/2-1,1);
mbnCtr=.71;
mbnDZ=.15;
mbnR=.38;
ctr=ceil(n/2+1);
mbn=fuzzymask(n,3,[mbnR mbnR mbnDZ]*n,n/20,[ctr ctr mbnCtr*n+1]);
%ShowSections(ms+5*max(mbn,blurMsk),[ctr,ctr,round(mbnCtr*n+1)],45);

% With membrane disc
mbnDZ=.26;
mbnDisc=fuzzymask(n,2,mbnR*n,n/20);
zs=round(n*(mbnCtr-mbnDZ/2)):round(n*(mbnCtr+mbnDZ/2));
mbn=zeros(n,n,n,'single');
for i=zs
    mbn(:,:,i)=mbnDisc;
end;
mbnMask=min(1,max(0,GaussFilt(max(mbn,blurMsk),.03)));
ShowSections(ms+5*mbnMask.*sphere,[ctr,ctr,round(mbnCtr*n+1)],45);
if numel(discMaskName)>0
    WriteMRC(mbnMask,pixA,[outPath discMaskName]);
    disp([discMaskName ' written.']);
end;


