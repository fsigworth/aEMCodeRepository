% reEMClassify2.m

% Parameters

makeFakeData=0;
nIters=4;
flagSaveFigures=false;
flagWriteClasses=true;
removeRings=1;

% inPath='/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1/Stack/';
% cd(inPath);
inPath='';

% outPath='/Users/fred/EMWork/Simulations/Kv/Reconstructions/';
outPath='';

useParFor=1;
nSlices=4;
doShowParticles=1;
mapName='KvMap.mat';  % for making fake data

% Working image sizes
% n=64;
nDs=32;    % downsampled box size
nCrop=24;  % final working size

disMag=1;  % display magnifications
disMag2=2;

% nDs=64;
% nCrop=64;  % final working size
% disMag=1;

% Define angle range

nGammas=10;
angleSteps=[5 5 360/nGammas]; % step sizes in degrees
angleLimits=[-10 20 360];  % starting alpha, beta; ending gamma
isos=[0 1];

doNormalize=0;
maxShift=3; % Maximum translation searched
nt=2*maxShift+1;

useEM1=false;  % Use the FFT-based EM step instead of the fast EM2 step

s0=.003;  % reference amplitude factor

sigmaC=3;
sigmaG=3;

imgAmp=1;
sigmaN=1;

nSimMicrographs=1000;
pSimIso=.8;
simImgScale=.002;

rng(1);  % initialize the random-number generator

if makeFakeData
    disp('Making fake data');
    [si0,imgs0]=reMakeFakeData(n,mapName,angleSteps,angleLimits,pSimIso,...
        simImgScale,nSimMicrographs);
else
    % ----------------Get the stack files---------------
    [si0, imgs0, names]=reLoadStackFiles;
    for i=numel(names)
        disp(names{i});
    end;
end;

[si0,imgs0]=rsStackSplit(si0.activeFlags(:,end),si0,imgs0);

% af=false(8006,3);
% nImgs=numel(si0.miIndex)
% load('aFlags1.mat')
% af(si0.origParticle,1)=true;
% load('aFlags2.mat')
% [si0,imgs0]=rsStackSplit(aFlags2,si0,imgs0);
% af(si0.origParticle,2)=true;
% 
% load('aFlags3.mat')
% [si0,imgs0]=rsStackSplit(aFlags3,si0,imgs0);
% af(si0.origParticle,3)=true;
% return

[si, imgs]=rsStackDownsample(si0,imgs0,nDs);
[si, imgs]=rsStackCrop(si,imgs,nCrop);

n0=size(imgs,1)
nImgs=size(imgs,3)

rotMask=fuzzymask(n0,2,n0/2,n0/20);

%%
%  We now have imgs and the si struct initialized.
if doShowParticles
    disp('Showing particles');
    figure(1);
    ImagicDisplay3(imgs);
    drawnow;
end;
%%
% ----------Set up the Run Info structure-----------
ri=struct;
ri.radius=floor(n0/2)-maxShift;
ri.nTrans=nt;
ri.softMask=fuzzymask(n0,2,ri.radius*.9,ri.radius*.2);
% Pick nominal angles for the references
[ri, refAngles]=reSetRefAngles(angleSteps,angleLimits,isos,false,ri);

nRefs=size(refAngles,1);
disp(['nAngles: ' num2str(ri.angleN) '  nRefs: ' num2str(nRefs)]);
disp(['shift range: ' sprintf('%d x %d',nt,nt)]);

% nAlphas=numel(ri.alphas);
nAlphasI=numel(ri.alphasI);
nIsos=numel(ri.isos);
% nRefs=ri.angleN(2)*ri.angleN(3);
nTA=nt^2*nAlphasI;

nRefs
refs=GaussFilt(randn(n0,n0,nRefs),.1,nRefs>1);  % random starting refs

% Set up the model info
moi.sigmaN=sigmaN;
moi.sigmaC=sigmaC;
moi.sigmaG=sigmaG;
moi.b0=0;
moi.pVols=1;
moi.imgAmps=imgAmp*ones(nImgs,1);
moi.a=mean(moi.imgAmps);
moi.pRefs=1/nRefs*ones(nImgs,nRefs,1,'single');


%%

tic
uImgs=zeros([n0 n0 nAlphasI nImgs],'single');

%---------------parfor loop for rotating images----------------
if useParFor
    disp('Rotating images -parallel');
    parfor i=1:nImgs
        uImgs(:,:,:,i)=rsRotateImage(imgs(:,:,i).*rotMask,-ri.alphasI);
    end;
else
    disp('Rotating images');
    for i=1:nImgs
        uImgs(:,:,:,i)=rsRotateImage(imgs(:,:,i).*rotMask,-ri.alphasI);
    end;
end;
toc


%% set up slicing the dataset
nST=nSlices;
nSliceImgs=ceil(nImgs/nST);

ro1=cell(nST,1);
imgBestMatch1=cell(nST,1);
classMeans1=cell(nST,1);
classNorms1=cell(nST,1);
si1=cell(nST,1);
rawLogProbs1=cell(nST,1);

% %%
sliceFlags=false(nImgs,nST);
    
    for is=1:nST
        sliceFlags(is:nST:nImgs,is)=true;
    end;

% refs=imgs;
% imgs=imgs+0.1*randn(size(imgs));

%
iter=0;
%%
while iter<nIters
    %%
    iter=iter+1

    %     Normalize the ref amplitudes to unity variance
    refsNorm=refs;
    rSigma=sqrt(refs(:)'*refs(:)/numel(refs));
    refsNorm=refsNorm/rSigma;

    % ------------------------ Do the EM step ---------------------------
    % ------------------------------------------------------------------
%                 ,rawMeans1,rawNorms1,rawLogProbs1{is}...
    pause(0.1);
    tic
    if useParFor
            disp(['parallel EM: slices=' num2str(nST)]);
        parfor is=1:nST
            [si1{is},uImg1]=rsStackSplit(sliceFlags(:,is),si,uImgs);
            mo1=reModelSplit(moi,si1{is}.pastParticle);
            [classMeans1{is}, classNorms1{is},...
                ro1{is},imgBestMatch1{is}]=reEMStep22(uImg1,refs,si1{is},ri,mo1);
        end;        
    else
        for is=1:nST
            [si1{is},uImg1]=rsStackSplit(sliceFlags(:,is),si,uImgs);
            mo1=reModelSplit(moi,si1{is}.pastParticle);
                [classMeans1{is}, classNorms1{is},...
                ro1{is}, imgBestMatch1{is}...
            ]=reEMStep22(uImg1,refs,si1{is},ri,mo1);
        end;
    end;
    toc
    % ------------------------------------------------------------------
    % ------------------------------------------------------------------
            
    % -----Accumulate slices----------
    
    
    roi=struct;
    roi.imgAmps=zeros(1,nImgs,'single');
    roi.pTransUnrot=zeros(nt,nt,nImgs,'single');
    roi.pTrans=zeros(nt,nt,nImgs,'single');
    roi.varNs=zeros(1,nImgs,'single');
    roi.logPX=zeros(1,nImgs,'single');
    roi.pRefs=zeros(nRefs,1,nImgs,'single');
    roi.pAlphas=zeros(nAlphasI,nImgs,'single');
    roi.imgClass=zeros(1,nImgs,'single');
    roi.imgTA=zeros(3,nImgs,'single');
    roi.pVols=zeros(1,nImgs,'single');
    imgBestMatch=zeros(n0,n0,nImgs,'single');
    
%     rawLogProbs=struct;
%     rawLogProbs.priors=zeros(nTA,nRefs,nImgs,'single');
%     rawLogProbs.pXs=zeros(nTA,nRefs,nImgs,'single');
    
    for is=1:nST
        inds=si1{is}.pastParticle;
        %  Merge the accumulated parameters
        roi.imgAmps(1,inds)=ro1{is}.imgAmps;
        roi.pTransUnrot(:,:,inds)=ro1{is}.pTransUnrot;
        roi.pTrans(:,:,inds)=ro1{is}.pTrans;
        roi.varNs(1,inds)=ro1{is}.varNs;
        roi.logPX(1,inds)=ro1{is}.logPX;
        roi.pRefs(:,:,inds)=ro1{is}.pRefs;
        roi.pAlphas(:,inds)=ro1{is}.pAlphas;
        roi.imgClass(1,inds)=ro1{is}.imgClass;
        roi.imgTA(:,inds)=ro1{is}.imgTA;
        roi.pVols(:,inds)=ro1{is}.pVols;
        imgBestMatch(:,:,inds)=imgBestMatch1{is};
        
%         rawLogProbs.priors(:,:,inds)=rawLogProbs1{is}.priors;
%         rawLogProbs.pXs(:,:,inds)=rawLogProbs1{is}.pXs;
        
    end;
    
    
    pause(0.1);
    
    % accumulate the fourier projections

        for is=1:nSlices;
            if is==1
                classMeans=classMeans1{1};
                classNorms=classNorms1{1};
            else
                classMeans=classMeans+classMeans1{is};
                classNorms=classNorms+classNorms1{is};
            end;
        end;
        classNorms=circshift(classNorms,ceil([n0/2 n0/2 0]));  %2D ifftshift
        
        reShowLatentVars(imgs,refs,ri,roi,iter,[4 10],'');
            drawnow;
    % ----------------Do the Fourier reconstruction----------------
    
    disp('Reconstructions');
    k=1;
    newRefs=real(ifft2(fft2(classMeans)./(classNorms+k)));
    
    pRefMin=1e-4;
    
    % Assign new model values
    moi.sigmaN=sqrt(mean(roi.varNs));
    moi.a=median(roi.imgAmps);
%     moi.imgAmps=0*roi.imgAmps+median(roi.imgAmps);
    moi.imgAmps=roi.imgAmps;
    moi.pRefs=max(shiftdim(roi.pRefs,2),pRefMin);
    
    oldRefs=refs;
    refs=newRefs;
    
    removeRings=iter>2;
    if removeRings
        disp('Remove rings');
        mRefs=reMakeMatchedRefs(refs,si,ri,roi);
        %%
        radii=[-30 -20 -10 0 10 20 30]/si.pixA;
        widths=radii*0+max(1,14/si.pixA);
        fitMask=(0.8*ri.softMask+0.2);
        fits=reFitVesicleRings(si,imgs-mRefs,fitMask,radii,widths);
        corrImgs=imgs-fits;
        %%
        disp('Rotating images');
        
        alphas=ri.alphasI;
        parfor i=1:nImgs
            uImgs(:,:,:,i)=rsRotateImage(corrImgs(:,:,i).*rotMask,-alphas);
        end;
    end;

    
    moi
    %
    if flagWriteClasses
            outName=sprintf('%sCls%di.mrc',outPath,iter);
            disp(outName);
    end;
    
    figure(10);
    %     subplot(nr,nc,3);
    %     title(['Iteration ' num2str(iter)]);
    if flagSaveFigures
        set(gcf,'paperpositionmode','auto');
        outName=sprintf('%sImgsIter%d.jpg',outPath,iter);
        print('-djpeg','-r200',outName);
        save([outPath 'fsc.mat'], 'fsc');
    end;
    
    drawnow;
    
        figure(1);
        ImagicDisplay3(newRefs);
        drawnow;
    %     q=reshape(roi.pRefs,nRefs*nVols,nImgs);
    % [vals, roi.classInds]=max(q);
    % roi.classInds=roi.classInds(:);
end; % while iter
%%

pTr=reshape(roi.pTrans,nt^2,nImgs);
[~,mxi]=max(pTr);
[tx,ty]=ind2sub([nt nt],mxi);
tx=tx'-maxShift-1;
ty=ty'-maxShift-1;
[~,mxa]=max(roi.pAlphas);
alphasRot=ri.alphasI(mxa)';
[~,mxr]=max(roi.pRefs,[],1);
mxr=mxr(:);
cls=cell(nRefs,1);
for i=1:nRefs
    cls{i}=single(find(mxr==i));
end;

% Make the rotated and aligned stack
% rImgs=zeros(n0,n0,nImgs,'single');
nCls=zeros(nRefs,1);
clsLists=cell(nRefs,1);
clsMeans=zeros(n0,n0,nRefs,'single');
for i=1:nImgs
%     urotInd=find(ri.alphas==-alphasRot(i));
%     if numel(urotInd)~=1
%         error('couldn''t match the angle');
%     end;
      urotInd=mxa(i);
    rImg=circshift(uImgs(:,:,urotInd,i),-[tx(i) ty(i)]);
%     imacs(rImg);
%     title(i);
%     pause(0);
    j=mxr(i);
    clsMeans(:,:,j)=clsMeans(:,:,j)+rImg;
    nCls(j)=nCls(j)+1;
end;
for j=1:nRefs
    if nCls(j)>0
        clsMeans(:,:,j)=clsMeans(:,:,j)/sqrt(nCls(j));
    end;
end;
figure(3);
ImagicDisplay3(ExpandImage(clsMeans,disMag2));


partCls=mxr;
        clsName=[names{end}(1:end-7) 'n' num2str(n0) 'c' num2str(nRefs) 'cls.mat'];

save(clsName,'partCls','clsLists','clsMeans');


return

%% Test code
figure(1);
% ImagicDisplay2(imgs,disMag);
ImagicDisplay3(newRefs);
imIndex=1;
qp0=rawLogProbs.priors(:,imIndex,1);
qp=reshape(qp0,nt,nt,nAlphas*nIsos);
qpn=reshape(rawLogProbs.priors(:,:,1),nt,nt,nAlphas*nIsos*nRefs);
qp1=qpn(6,:,4:9:nAlphas*nIsos*nRefs);
qx0=rawLogProbs.pXs(:,imIndex,1);
qx=reshape(qx0,nt,nt,nAlphas,nIsos);


figure(8);
SetGrayscale;
imovie(qx,.1);
qx1=qx(:,:,end);

figure(9);
SetGrayscale;

% imovie(qp,.1);
plot(sect(rawLogProbs.priors(:,:,1)));

figure(11);
SetGrayscale;
imacs(rawLogProbs.priors(:,:,1));
plot(sect(rawLogProbs.priors(:,:,1)));
% imacs(rawLogProbs.pXs(:,:,1));