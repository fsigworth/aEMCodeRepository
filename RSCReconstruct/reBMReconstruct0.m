% reBMReconstruct0
% Best-match assignment of angles based on the roi structure from
% reEMReconstruct.  
stackSuffix='ustack.mrc';
useAltImgs=1;  % reconstruct with unsub images
n=0; % use original value
disp('Loading roi and ri');
% load('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon112j4/i40a_roi.mat')
% load('/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/Reconstructions/Recon112j4/ri.mat')
% iVol=1;
% baseDir='/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/';
%reconDir='Reconstructions/Recon72k1/';
% cd(baseDir);
% load([reconDir 'i28a_roi.mat'])
% load([reconDir 'ri.mat'])
% iVol=4;
% newSiName='sq03_1v2fp256tsi.mat';

% reconDir='Reconstructions/Recon112j4/';
% cd(baseDir);
% load([reconDir 'i40a_roi.mat'])
% load([reconDir 'ri.mat'])
% iVol=1;  % not so good looking, but has 88% of particles
% ri0=ri;  % keep the original

% baseDir='/Users/fred/EMWork/Hideki/151117/KvLipo80slot3/';
% reconDir='Reconstructions/Recon72k4v4x/';
% cd(baseDir);
% load([reconDir 'i48a_roi.mat'])
% load([reconDir 'ri.mat'])
% iVol=4;  % 3 is biggest class; 4 maybe looks nicer
% ri0=ri;  % keep the original

% rootDir='/Users/fred/EMWork/Hideki/';
rootDir='/Volumes/D213/EMWork/';
baseDir=[rootDir '151117/KvLipo80slot3/'];
reconDir='Reconstructions/Recon112s10v2r1/';
cd(baseDir);
load([reconDir 'i29a_roi.mat'])
load([reconDir 'riExpandedMask.mat'])
iVol=2;  %
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

rip=ri;  % subtracted particles

rip.nGroups=1;
if n>0 % we want to re-scale
    rip.nCurrent=n;  % downsampling
    ds=ri.nCurrent/n;
    rip.softMask=DownsampleGeneral(rip.softMask,n);
else
    ds=1;  % no rescaling
end;
iTwin=1;
iGroup=1;

%%
disp('Loading and scaling the particle images');

[pSi,pImgs,pzImgs]=reLoadStackGroup(rip,si,iTwin,iGroup,useAltImgs);
if useAltImgs
    pImgs=pzImgs;
end;
%%
nim=size(pImgs,3);
disp('Rotating images');
for i=1:nim
    pImgs(:,:,i)=pImgs(:,:,i).*rip.softMask;
end;
tic
pRotImgs=rsRotateImage(pImgs,-rip.alphas(bmAIs));
toc
%


% Fvs=AccumAndInsert(rotImgs,pVs,
% gRotImgs=rotVesImgs;
isos=[0 1];
n1=size(pRotImgs,1);

pTransRotImgs=zeros(n1,n1,nim,'single');

% do translations
for i=1:nim
        imgAmp=roi.imgAmps(i);
        pP=FourierShift(rip.nCurrent,-[tx(i),ty(i)]/ds+(nt/2+1));  % not downsampled
        pTransRotImgs(:,:,i)=real(ifftn(fftn(pRotImgs(:,:,i)).*pP));
end;

disp('Accumulating classes');
isVol = (pVs(iVol,:)>.5);
isIsos=false(1,nim);
for i=1:numel(isos)
    isIsos=(bmIsos==isos(i)) | isIsos;
end;

isActiveImg=isVol & isIsos;
imgAmps=roi.imgAmps;


% Save all the data
disp('saving...');
save([reconDir 'ParticleStack112v2U.mat'],'pTransRotImgs');
save([reconDir 'ParticleVars112v2U.mat'], 'rip','pSi','bmRefs','isActiveImg','imgAmps');
% rip: run info structure for subtracted particles
% pSi: stack info structure
% bmRefs: index of projection directions (into rip.angles) for image
% isActiveImg: booleans of active images






%%-------------------------------
% clear;
% disp('loading...');
% load ParticleStack.mat
% load Recon72v4.mat

[n1,~,nim]=size(pTransRotImgs);
nRefs=size(rip.angles,1);
activeImgs=find(isActiveImg);
pClassMeans=zeros(n1,n1,nRefs,'single');
pClassNorms=zeros(n1,n1,nRefs,'single');

for i=activeImgs
        imgAmp=imgAmps(i);
        irv=bmRefs(i);  % best-match reference index        
        pCtf=ifftshift(pSi.ctfs(:,:,pSi.miIndex(i)));
        pCImg=real(ifftn(pCtf.*fftn(pTransRotImgs(:,:,i)))); % CTF-filtered image
        pClassMeans(:,:,irv)=pClassMeans(:,:,irv)+pCImg*imgAmp;  % real space img sum
        pClassNorms(:,:,irv)=pClassNorms(:,:,irv)+pCtf.^2*imgAmp^2; % Fourier, zero-centered
end;
disp([num2str(numel(activeImgs)) ' images accumulated from ' num2str(nim)]);
%
useParFor=1;

pRealNorms=zeros(n1,n1,nRefs,'single');
pProjs=zeros(n1,n1,nRefs,2,'single');
for i=1:nRefs
    pRealNorms(:,:,i)=fftshift(real(ifftn(pClassNorms(:,:,i))));
end;
pProjs(:,:,:,1)=pClassMeans;
pProjs(:,:,:,2)=pRealNorms;
disp('Fourier insertion--particles');
tic
pFvs=reFourierInsertion2(pProjs,rip.angles,rip.symmetry,useParFor);
toc
%%         [moi,vols]=reGroupReconstruct(ri,si,iter,iTwin,roi,oldMoi,gFvs,logs);
%
disp('Volume normalization');
pK=mean(imgAmps.^2)*rip.kFactor;  % mult by fsc/(1-fsc) says Sjors
pK=mean(roi.varNs)*.1;
pVol=single(rsNormalizeReconstruction(pFvs(1),pFvs(2),pK));
%
figure(10);
ShowSections(pVol,[],45);
%%
angles=ri.angles;
save([reconDir 'ReconV2U.mat'], 'pVol','pFvs','pProjs','angles', 'pK');

