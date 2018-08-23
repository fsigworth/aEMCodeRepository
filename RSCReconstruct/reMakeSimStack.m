% reMakeSimStack.m

% mapName='3KG2mapsub2.9A.mat';
mapName='KvMap.mat';
gammaStepFactor=2;

outPath='/Users/fred/EMWork/Hideki/160622Sim/Stack/';
outPath='/fastscratch/fjs2/1606225Sim/Stack';

outPostfix='';

% Working image sizes
n=128;  % final image size
nc=NextNiceNumber(1.2*n);  % padded image for CTF etc.
pixA=4.988/2;
symmetry=4;
nMicImgs=40;
nMicrographs=200;
sigmaN=.7;
imgAmp=.001;

B0=40;
B1=0;
ctfAlpha=.02;

betaMin=20;
defs=[2 3.5];  % Defocus range
vesRs=[350 500]/pixA; % Vesicle radius range
pIso=.5;

lambda=.025;

mi0=meCreateMicrographInfoStruct14;
mi0.pixA=1.25;
mi0.camera=5;
mi=mi0;
mi(nMicrographs).pixA=mi0.pixA;

% Angle and shift search

sigmaC=0;
sigmaGr=0;
sigmaGb=0;
biasBob=0;
biasAlphaBeta=[0 0];
biasClick=[0 0];

nVols=1;  % only implemented 1 volume at present

[vols, mbnOffsetA]=arGetRefVolumes(pixA,n,mapName,nVols);
mbnOffset=mbnOffsetA/pixA;

nImgs=nMicImgs*nMicrographs;
nImgs

si=struct;
si.miIndex=uint16(1+floor(0:1/nMicImgs:nMicrographs-1/nMicImgs)');
si.miParticle=repmat((1:nMicImgs)',nMicrographs,1);
si.alpha0=single(360*rand(nImgs,1)-180);
si.rVesicle=single(vesRs(1)+(vesRs(2)-vesRs(1))*rand(nImgs,1));
si.pixA=pixA;
si.mbnOffset=mbnOffsetA/pixA;
si.weights=[1 1];
si.ctfs=zeros(n,n,nMicrographs,'single');
si.sVesicle=0.003*ones(nImgs,1,'single');
si.mi=mi;
si.activeFlags=true(nImgs,1);
si.activeFlagSet={['reMakeSimStack  ' date]};

sm.defoci=defs(1)+(defs(2)-defs(1))*rand(nMicrographs,1);

ctfsc=zeros(nc,nc,nMicrographs,'single');
for j=1:nMicrographs
    d=sm.defoci(j);
    si.ctfs(:,:,j)=abs(CTF(n,pixA,lambda,d,2,B0+B1*d,ctfAlpha));
    ctfsc(:,:,j)=abs(CTF(nc,pixA,lambda,d,2,B0+B1*d,ctfAlpha));
end;
sm.isos=rand(nImgs,1)<pIso;
doIsoFlip=0;
% [ri, sm.angles]=reSetRefAngles(angleSteps,angleLimits,doIsoFlip,false,ri);

sm.angles=rsSphereAngles(nImgs,symmetry,betaMin,gammaStepFactor);
sm.rocks=sigmaGr*randn(nImgs,2)+repmat(biasAlphaBeta,nImgs,1);  % alpha,beta rocking
sm.bobs=sigmaGb*randn(nImgs,1)+biasBob;
sm.clicks=sigmaC*randn(nImgs,2)+repmat(biasClick,nImgs,1);
sm.mbnOffset=mbnOffset;

% Make the images -------------
[imgs0,si.yClick,sm.partAngles]=rsMakeFakeImages2(vols,nc,si.rVesicle,sm);
msk=fuzzymask(nc,2,0.45*n,0.1*n);
imgs1=zeros(nc,nc,nImgs,'single');
for i=1:nImgs
    c=ctfsc(:,:,si.miIndex(i));
    imgs1(:,:,i)=msk.*real(ifftn(fftn(imgs0(:,:,i)).*ifftshift(c)));
end;
%%
imgsC=imgAmp*Crop(imgs1,n,1);
imgs=imgsC+sigmaN*randn(n,n,nImgs);
%
figure(1);
imagsar(BinImage(imgs,4));

outName=sprintf('Simn%dsiz%02dsigma%02dclk%02d%s',nImgs,n,round(sigmaN*100),round(sigmaC),outPostfix);
disp(['Writing: ' outName]);
siName=[outPath outName 'tsi.mat'];
disp(siName);
save(siName,'si');
smName=[outPath outName 'tsm.mat'];
disp(smName);
save(smName,'sm');

stackName=[outPath outName 'tstack.mrc'];
disp(stackName);
WriteMRC(imgs,pixA,stackName);

% Make the low-noise stack
sigmaN0=sigmaN*.02;
outName0=sprintf('Simn%dsiz%02dsigma%02dclk%02d%s',nImgs,n,round(sigmaN0*100),round(sigmaC),outPostfix);
siName0=[outPath outName0 'tsi.mat'];
disp(siName0);
save(siName0,'si');  % identical file to the other one.

imgs0=imgsC+sigmaN0*randn(n,n,nImgs);
stackName0=[outPath outName0 'tstack.mrc'];  % no noise
disp(stackName0);
WriteMRC(imgs0,pixA,stackName0);
