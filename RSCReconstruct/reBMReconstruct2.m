% reBMReconstruct2
% Best-match reconstruction based on the roi structure from
% reEMReconstruct.  in Version 2 it is used to reconstruct in a larger volume.
useAltImgs=1;
stackSuffix='ustack.mrc';
ds=2;  % downsampling of vesicle images

disp('Loading roi and ri');
% load('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon112j4/i40a_roi.mat')
% load('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon112j4/ri.mat')
% iVol=1;
baseDir='/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/';
reconDir='Reconstructions/Recon72k1/';
cd(baseDir);
load([reconDir 'i28a_roi.mat'])
load([reconDir 'ri.mat'])
iVol=4;
newSiName='sq03_1v2fp256tsi.mat';
% 
% load('/Users/fred/EMWork/Hideki/140625n/Reconstructions/Recon112j6/i40a_roi.mat')
% load('/Users/fred/EMWork/Hideki/140625n/Reconstructions/Recon112j6/ri.mat')
% cd('/Users/fred/EMWork/Hideki/140625n')
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

% Modify ri to point to the new name
ri.siName=newSiName;
disp(['Loading the new si file ' ri.siName]);
load([ri.stackPath ri.siName]);  % get si
origN=size(si.ctfs,1);

us=origN/origN0  % size increase
ds
%%
ri.nCropU=us*ri0.nCropU;
ri.nFinal=us*ri0.nFinal/ds;
n2=ri.nFinal;
ri.nSequence(:,1)=ri0.nSequence(:,1)*us/ds;
ri.nCurrent=ri0.nCurrent*us/ds;
ri.volMask=Downsample(ri0.volMask,ri.nCurrent);  % scale up
ri.softMask=Downsample(ri0.softMask,ri.nCurrent); % particle mask
ri.fscMask=Downsample(ri0.fscMask,ri.nCurrent);

ri.nGroups=1;
iTwin=1;
iGroup=1;

disp('Loading and scaling the images');
% stackSet=[0 1];  % just the alt stack
stackSet=[1 1];
[tSi,tImgs,tAltImgs,groupSize,twinInds]=reLoadStackGroup(ri,si,iTwin,iGroup,stackSet);
tImgs=tAltImgs-tImgs;
clear tAltImgs

%%
nim=size(tImgs,3);

disp('Rotating the images');
tic
for i=1:nim
    tImgs(:,:,i)=tImgs(:,:,i).*ri.softMask;
end;
rotVesImgs=rsRotateImage(tImgs,-ri.alphas(bmAIs));
% gRotImgs0=rsRotateImage(gImgs,-ri.alphas(bmAIs));
toc
%%
gRotImgs=rotVesImgs;
isos=[0 1];

% n=size(gRotImgs,1);
classMeans=zeros(n2,n2,nRefs,'single');
classNorms=zeros(n2,n2,nRefs,'single');
transRotImgs=zeros(n2,n2,nim,'single');
disp('Accumulating classes');
nacc=0;
for i=1:nim
    if pVs(iVol,i)>.5 && any(bmIsos(i)==isos)% this volume, correct iso
        P=FourierShift(ri.nCurrent,-[tx(i),ty(i)]+(nt/2+1));
%         P=1;
        ctf=ifftshift(tSi.ctfs(:,:,tSi.miIndex(i)));
        transRotImgs(:,:,i)=real(ifftn(ctf.*fftn(gRotImgs(:,:,i)).*P));  % ifftshift?
        % pixMeans=eImg_tai_rv*imgAmp
        imgAmp=roi.imgAmps(i);
        irv=bmRefs(i);
        classMeans(:,:,irv)=classMeans(:,:,irv)+transRotImgs(:,:,i)*imgAmp;  % real space
        classNorms(:,:,irv)=classNorms(:,:,irv)+ctf.^2*imgAmp^2; % Fourier, zero-centered
        nacc=nacc+1;
    end;
end;
disp([num2str(nacc) ' images accumulated from ' num2str(nim)]);
%
% prune references
weights=sum(reshape(classNorms,n2^2,nRefs));
% % thresh=Percentile(weights,.5);
thresh=0;
activeRefs=weights(:)>thresh;

angles=ri.angles(activeRefs,:);
classMeans=classMeans(:,:,activeRefs);
classNorms=classNorms(:,:,activeRefs);
naRefs=sum(activeRefs);
disp([num2str(naRefs) ' active refs out of ' num2str(nRefs)]);

%%
useParFor=1;
    realNorms=zeros(n2,n2,naRefs,'single');
    projs=zeros(n2,n2,naRefs,2,'single');
    
    for i=1:naRefs
        realNorms(:,:,i)=fftshift(real(ifftn(classNorms(:,:,i))));
    end;
    projs(:,:,:,1)=classMeans;
    projs(:,:,:,2)=realNorms;
disp('Fourier insertion');
    Fvs=reFourierInsertion(projs,angles,ri.symmetry,useParFor);
% %         [moi,vols]=reGroupReconstruct(ri,si,iter,iTwin,roi,oldMoi,gFvs,logs);
disp('Volume normalization');

k=mean(roi.imgAmps.^2)*ri.kFactor;  % mult by fsc/(1-fsc) says Sjors
%%
% Accumulate the Fourier volumes
disp(['k ' num2str(k)]);
    v=single(rsNormalizeReconstruction(Fvs(1),Fvs(2),k));
moi.refVols=reNormalizeModels(ri,v);

figure(10);
ShowSections(Downsample(v,size(v)/2),[],45);
save([reconDir 'vIntermedSym4'], 'v', 'Fvs', 'classMeans', 'projs', 'angles');
return
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

