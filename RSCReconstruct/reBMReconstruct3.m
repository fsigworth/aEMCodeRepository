% reBMReconstruct3
% Best-match reconstruction based on the roi structure from
% reEMReconstruct.  in Version 3 it is used to reconstruct particles and
% membranes separately
useAltImgs=1;  % use membrane images
stackSuffix='ustack.mrc';
dsv=2;  % downsampling of vesicle images

disp('Loading roi and ri');

% 151117 dataset
% baseDir='/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/';
% reconDir='Reconstructions/Recon112k6v0/';
baseDir='/Volumes/D213/EMWork/151117/KvLipo80slot3/';
reconDir='Reconstructions/Recon72s10v4r1/';
cd(baseDir);
load([reconDir 'i28a_roi.mat'])
load([reconDir 'ri.mat'])
iVol=4;  %
newSiName='sq03_1v2fp256tsi.mat';
%
% baseDir='/Users/fred/EMWork/Hideki/140625/';
% reconDir='Reconstructions/Recon112jl_mms_copy/';
% cd(baseDir);
% load([reconDir 'i28a_roi.mat'])
% load([reconDir 'ri.mat'])
% newSiName='sq10_350v2fp256tsi.mat';

ri0=ri;  % keep the original


% Pick up the best-match angles and translations
disp('Getting best-match parameters');
[nAlphas,nIOs,~]=size(roi.pAlphas);
[nGammas,nBetas,nVols,nim]=size(roi.pRefs);
nRefs=nGammas*nBetas;
nimAF=sum(ri.activeFlags);
nt=sqrt(size(roi.pTrans,1));

% roi.pRefs[gammas, betas, vols, imgs]
pRVs=reshape(roi.pRefs,nRefs,nVols,nim);
% Get the volume probability, for particle selection
pVs=shiftdim(sum(sum(roi.pRefs,2),1),2); % result is nvols x nim
% Get the best-matching ref index
pRs=squeeze(pRVs(:,iVol,:));  % prob(refs,imgs) for selected volume
[~,bmRefs]=max(pRs,[],1); % best ref for each image

% Best-matching alpha and iso
pAIs=reshape(roi.pAlphas,nAlphas*nIOs,nim);
[~,bmAIs]=max(pAIs,[],1);

bmIsos=bmAIs>nAlphas;

inBoundImgs=roi.inBoundImgs;

% Best-matching translation
pTrans2=reshape(roi.pTrans,nt,nt,nim);
tx=zeros(nim,1);
ty=zeros(nim,1);
for i=1:nim
    [~,tx(i),ty(i)]=max2di(pTrans2(:,:,i));
end;

% Load the si files
disp(['Loading the old si file ' ri0.siName]);
load([ri0.stackPath ri0.siName]);
si0=si;  % keep the original stack
origN0=size(si0.ctfs,1);

% Make a new ri file for the big vesicle volume.
rip=ri;  % subtracted particles
rip.siName=newSiName;
rip.nGroups=1;
n1=rip.nCurrent;

riv=ri;  % vesicles
riv.siName=newSiName;  % point to the new volume
disp(['Loading the new si file ' riv.siName]);
load([riv.stackPath riv.siName]);  % get si
origN=size(si.ctfs,1);
si.activeFlags=si0.activeFlags;  % make these consistent.

usv=origN/origN0  % size increase of new relative to old
dsv


%%  Modify the riv structure
%   size increase by us, but downsampling by ds
riv.nCropU=usv*ri0.nCropU;
riv.nFinal=usv*ri0.nFinal/dsv;
n2=riv.nCurrent;
riv.pixA=ri0.pixA*dsv;
riv.nSequence(:,1)=ri0.nSequence(:,1)*usv/dsv;
riv.nCurrent=ri0.nCurrent*usv/dsv;
riv.volMask=Downsample(ri0.volMask,riv.nCurrent);  % scale up
riv.softMask=Downsample(ri0.softMask,riv.nCurrent); % particle mask
riv.fscMask=Downsample(ri0.fscMask,riv.nCurrent);
riv.nGroups=1;
iTwin=1;
iGroup=1;

disp('Loading and scaling the vesicle images');
stackSet=[1 1];
[vSi,tImgs,tAltImgs]=reLoadStackGroup(riv,si,iTwin,iGroup,stackSet);
%
disp('Loading and scaling the particle images');
vImgs=tAltImgs-tImgs;  % vesicles alone
% vImgs=tAltImgs;  %%%%%%%%% vesicles and particles
stackSet=[1 0];
[pSi,pImgs,pzImgs]=reLoadStackGroup(rip,si,iTwin,iGroup,stackSet);

%
nim=size(vImgs,3);
disp('Rotating images');
for i=1:nim
    vImgs(:,:,i)=vImgs(:,:,i).*riv.softMask;
    pImgs(:,:,i)=pImgs(:,:,i).*rip.softMask;
end;
%
tic
vRotImgs=rsRotateImage(vImgs,-riv.alphas(bmAIs));
toc
tic
pRotImgs=rsRotateImage(pImgs,-rip.alphas(bmAIs));
% gRotImgs0=rsRotateImage(gImgs,-ri.alphas(bmAIs));
toc
%     vRotImgs=vImgs; %%%%%%%
%     pRotImgs=pImgs;
%%
% Fvs=AccumAndInsert(rotImgs,pVs,
% gRotImgs=rotVesImgs;
volThreshold=0.7;
isos=[0 1];
% isos=0; % rso only
% n=size(gRotImgs,1);
vClassMeans=zeros(n2,n2,nRefs,'single');
vClassNorms=zeros(n2,n2,nRefs,'single');
vTransRotImgs=zeros(n2,n2,nim,'single');
pClassMeans=zeros(n1,n1,nRefs,'single');
pClassNorms=zeros(n1,n1,nRefs,'single');
pTransRotImgs=zeros(n1,n1,nim,'single');

disp('Accumulating classes');
nacc=0;
pAngles=rip.angles;
class=zeros(nim,1,'single');
for i=1:nim
    if pVs(iVol,i)>volThreshold && any(bmIsos(i)==isos) && inBoundImgs(i) % this volume, correct iso
        imgAmp=roi.imgAmps(i);
        irv=bmRefs(i);
        class(i)=irv;
        vP=FourierShift(riv.nCurrent,(-[tx(i),ty(i)]+(nt/2+1)) / dsv);
        vCtf=ifftshift(vSi.ctfs(:,:,vSi.miIndex(i)));
        vTransRotImgs(:,:,i)=real(ifftn(vCtf.*fftn(vRotImgs(:,:,i)).*vP));
        vClassMeans(:,:,irv)=vClassMeans(:,:,irv)+vTransRotImgs(:,:,i)*imgAmp;  % real space
        vClassNorms(:,:,irv)=vClassNorms(:,:,irv)+vCtf.^2*imgAmp^2; % Fourier, zero-centered
        
        pP=FourierShift(riv.nCurrent,-[tx(i),ty(i)]+(nt/2+1));  % not downsampled
        pCtf=ifftshift(pSi.ctfs(:,:,pSi.miIndex(i)));
        pTransRotImgs(:,:,i)=real(ifftn(pCtf.*fftn(pRotImgs(:,:,i)).*pP));
        pClassMeans(:,:,irv)=pClassMeans(:,:,irv)+pTransRotImgs(:,:,i)*imgAmp;  % real space
        pClassNorms(:,:,irv)=pClassNorms(:,:,irv)+pCtf.^2*imgAmp^2; % Fourier, zero-centered
        nacc=nacc+1;
    end;
end;
disp([num2str(nacc) ' images accumulated from ' num2str(nim)]);
%%
%
% prune references
weights=sum(reshape(vClassNorms,n2^2,nRefs));
% % thresh=Percentile(weights,.5);
thresh=0;
activeRefs=weights(:)>thresh;

angles=riv.angles(activeRefs,:);
vClassMeans=vClassMeans(:,:,activeRefs);
vClassNorms=vClassNorms(:,:,activeRefs);
pClassMeans=pClassMeans(:,:,activeRefs);
pClassNorms=pClassNorms(:,:,activeRefs);
naRefs=sum(activeRefs);
disp([num2str(naRefs) ' active refs out of ' num2str(nRefs)]);

%
useParFor=1;
vRealNorms=zeros(n2,n2,naRefs,'single');
vProjs=zeros(n2,n2,naRefs,2,'single');

for i=1:naRefs
    vRealNorms(:,:,i)=fftshift(real(ifftn(vClassNorms(:,:,i))));
end;
vProjs(:,:,:,1)=vClassMeans;
vProjs(:,:,:,2)=vRealNorms;
disp('Fourier insertion--vesicles');
tic
vFvs=reFourierInsertion2(vProjs,angles,riv.symmetry,useParFor);
toc

pRealNorms=zeros(n1,n1,naRefs,'single');
pProjs=zeros(n1,n1,naRefs,2,'single');
for i=1:naRefs
    pRealNorms(:,:,i)=fftshift(real(ifftn(pClassNorms(:,:,i))));
end;
pProjs(:,:,:,1)=pClassMeans;
pProjs(:,:,:,2)=pRealNorms;
disp('Fourier insertion--particles');
tic
pFvs=reFourierInsertion2(pProjs,angles,rip.symmetry,useParFor);
toc
%%         [moi,vols]=reGroupReconstruct(ri,si,iter,iTwin,roi,oldMoi,gFvs,logs);
%
disp('Volume normalization');

fc=.15
vK=mean(roi.imgAmps.^2)*riv.kFactor;  % mult by fsc/(1-fsc) says Sjors
disp(['vK ' num2str(vK)]);
f3n=(Radius3(pFvs(1).np1)/(pFvs(1).np*rip.pixA)/fc).^4;
H=(.01+f3n)./(1+f3n);
vKf=H.*vK;

vVol=single(rsNormalizeReconstruction(vFvs(1),vFvs(2),vKf));
vBigVol=Downsample(vVol,n2*dsv);

%

pK=vK;
fc=.15;
f3n=(Radius3(pFvs(1).np1)/(pFvs(1).np*rip.pixA)/fc).^4;
H=(.01+f3n)./(1+f3n);
pKf=H.*pK*1;

pVol=single(rsNormalizeReconstruction(pFvs(1),pFvs(2),pKf));
pBigVol=Crop(pVol,n1*usv);
tBigVol=vBigVol+0.6*pBigVol;
%
figure(10);
ShowSections(pVol,[],45);
% ShowSections(Crop(pBigVol+vBigVol*sqrt(8),n1),[],45);
%
save([reconDir 'P3ViorsoN1'], 'vVol', 'pVol', 'tBigVol','pFvs', 'vFvs','pProjs', 'vProjs','angles', 'pK');
save([reconDir 'P3ViorsoN1Imgs'],'pTransRotImgs', 'vTransRotImgs','class','pAngles','rip');
%%
load([reconDir 'P3ViorsoN1']);
% I find that Crop(pVol,144) has to be scaled down by 2 to match density of
% direct reconstruction of the entire 144^3 volume in reBMReconstruc.m

return
%% Make matching projections

% Make projections
refProjs=reMakeTemplates(pVol, angles);
%%
nc=9;
nProj=nc*ceil(size(angles,1)/nc);
mergedProjs=zeros(n1,n1,nc,nProj/nc*2);
rp=refProjs;
rp(1,1,nProj)=0;
rp=reshape(rp,n1,n1,nc,nProj/nc);
pc=pClassEsts(:,:,activeRefs);
pc(1,1,nProj)=0;
pc=reshape(pc,n1,n1,nc,nProj/nc);
mergedProjs(:,:,:,1:2:end)=rp;
mergedProjs(:,:,:,2:2:end)=pc;
mergedProjs=reshape(mergedProjs,n1,n1,nProj*2);
return
%%
% load([reconDir 'PUnsub']);
load([reconDir 'PAndV1']);
n1=size(pVol,1);
si.pixA=2.4;
tVol=1.7*Downsample(Crop(vVol,64),128)+Crop(pVol,128);
figure(2);
imagsar(tVol);
kMag=1.2;
q=load('/Users/fred/Structures/kv1.2/KvMap.mat');
kvMap=DownsampleGeneral(q.sim,n1,kMag*q.pixA/si.pixA);
rKvMap=ERotate3(kvMap,[0 pi/2 pi/2]);
figure(3);
imagsar(Crop(rKvMap,128));






%% Show class means
ncls=size(pClassMeans,3);
n=size(pClassMeans,1);
fc=.1;
nf=4;
fs=RadiusNorm(n)/rip.pixA;
pK=mean(roi.imgAmps.^2)*rip.kFactor;  % mult by fsc/(1-fsc) says Sjors
pKf=ifftshift(pK*(.001+(fs/fc).^nf)./(1+(fs/fc).^nf));
pClassEsts=pClassMeans*0;
for i=1:ncls
    pClassEsts(:,:,i)=real(ifftn(fftn(pClassMeans(:,:,i))./(pKf+pClassNorms(:,:,i))));
end;

figure(11);
imagsar(pClassEsts);








%%
cd(reconDir)
v0=ReadMRC('v0.mrc');
ImagicDisplay3(Crop(v0,96),0,1);

%%
figure(9);
text={'Right-side out particles';'Inside-out particles'};

active=~bmIsos;
for j=1:2
    inds=find(active);
    ni=numel(inds);
    rvals=zeros(ni,1);
    for i=1:ni
        ind=inds(i);
        rvals(i)=si.rVesicle(ind);
    end;
    subplot(2,1,j);
    hist(rvals*si.pixA,10:10:490);
    title(text{j});
    ylabel('Frequency');
    xlabel('Vesicle radius, Å');
    active=~active;
end;

