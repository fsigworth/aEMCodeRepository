% reEMClassify2.m

% Parameters

nIters=8;
flagSaveFigures=false;
flagWriteClasses=false;
%outPath='/Users/fred/EMWork/Simulations/Kv/Reconstructions/';
useParFor=1;
nSlices=16;
doShowParts=0;

% Working image sizes
nDs=64;
nCrop=64;  % final working size
% nDs=32;
% nCrop=32;  % final working size


% Define angle range
nGammas=5;
angleSteps=[2 10 360/nGammas]; % step sizes in degrees
angleLimits=[-10 10 360];  % starting alpha, beta; ending gamma
doIsoFlip=0;

% nGammas=1;
% angleSteps=[2 90 360/nGammas]; % no changes in beta
% angleLimits=[-12 90 360];  % 

doNormalize=1;
maxShift=8; % Maximum translation searched
nt=2*maxShift+1;

useEM1=false;  % Use the FFT-based EM step instead of the fast EM2 step

s0=.003;  % reference amplitude factor

sigmaC=20;
sigmaG=20;

imgAmp=1;
sigmaN=1;


    [fname, pa]=uigetfile('*si.mat','Select si files','multiselect','on');
    if isnumeric(pa)  % user clicked Cancel
        return
    end;
    if ~iscell(fname)
        fname={fname};
    end;
    cd(pa);
    outPath=AddSlash(pwd);
    
    
    %%
    % Assemble a complete stack
    allImgs=[];
    allSi=struct;
    for i=1:numel(fname)
        siName=fname{i};
        load(siName);
        p=regexp(siName,'si.mat');
        if numel(p)<1
            error(['No match found in name ' siName]);
        end;
        stackName=[siName(1:p(end)-1) 'stack.mrc'];
        disp(stackName);
        imgs=ReadEMFile(stackName);
        [allSi, allImgs]=rsStackConcatenate(si,imgs,allSi,allImgs);
    end;
    
    siName=fname{end};
    si=allSi;
    % Get rid of the original particle info
si.origParticle=uint32(1:numel(si.miIndex));

n0=size(allImgs,1);
nImgs=size(allImgs,3);

if doNormalize
    % Make an annulus for computing mean and variance
    %%%%% hardwired??
    ringo=0.48*n0;
    ringi=0.42*n0;
    msko=fuzzymask(n0,2,ringo,ringo/10);
    mski=fuzzymask(n0,2,ringi,ringi/10);
    annulus=msko-mski;
    
    % Compute mean and variance
    nann2=annulus(:)'*annulus(:);
    nann=sum(annulus(:));
    vars=zeros(nImgs,1);
    avgs=zeros(nImgs,1);
    for i=1:nImgs
        pix=allImgs(:,:,i).*annulus;
        avgs(i)=sum(pix(:))/nann;
        vars(i)=pix(:)'*pix(:)/nann2-avgs(i)^2;
    end;
    
    % Normalize variance to 1 in each image
    imgs=allImgs;
    for i=1:nImgs
        imgs(:,:,i)=(allImgs(:,:,i)-avgs(i))/sqrt(vars(i));
    end;
else
    imgs=allImgs;
end;

sigma=1;  % initialize the model parameter

[si, imgs]=rsStackDownsample(si,imgs,nDs);
[si, imgs]=rsStackCrop(si,imgs,nCrop);

%%

n0=size(imgs,1)
nImgs=size(imgs,3)


%  We now have imgs and the si struct initialized.
if doShowParts
    figure(1);
ImagicDisplay1(imgs);
drawnow;
end;

% Set up the run's 2d mask info
ri=struct;
ri.radius=floor(n0/2)-maxShift;
ri.nTrans=nt;
ri.softMask=fuzzymask(n0,2,ri.radius*.97,ri.radius*.1);
ri.symmetry=1;

ri=reSetRefAngles(angleSteps,angleLimits,doIsoFlip,1,ri);
nRefs=ri.angleN(2)*ri.angleN(3);

nRefs

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
uImgs=zeros([n0 n0 numel(ri.alphas) nImgs],'single');

%---------------parfor loop for rotating images----------------
if useParFor
    disp('Rotating images -parallel');
    parfor i=1:nImgs
        uImgs(:,:,:,i)=rsRotateImage(imgs(:,:,i),-ri.alphas);
    end;
else
    disp('Rotating images');
    for i=1:nImgs
        uImgs(:,:,:,i)=rsRotateImage(imgs(:,:,i),-ri.alphas);
    end;
end;
toc


%% iterate
nAlphas=numel(ri.alphas);

nST=nSlices;

nSliceImgs=ceil(nImgs/nST);

ro1=struct;
ro1.imgAmps=zeros(1,nSliceImgs,'single');
ro1.pTransUnrot=zeros(nt,nt,nSliceImgs,'single');
ro1.pTrans=zeros(nt,nt,nSliceImgs,'single');
ro1.varNs=zeros(1,nSliceImgs,'single');
ro1.pVols=ones(1,nSliceImgs,'single');
ro1.logPX=zeros(1,nSliceImgs,'single');
ro1.pRefs=zeros(nRefs,1,nSliceImgs,'single');
ro1.pAlphas=zeros(nAlphas,nSliceImgs,'single');
ro1.imgClass=zeros(1,nSliceImgs,'single');
ro1.imgTA=zeros(3,nSliceImgs,'single');

ro1(nST,1).imgAmps=[];  % extend the struct array.

imgBestMatch1=cell(nST,1);
classMeans1=cell(nST,1);
classNorms1=cell(nST,1);
si1=cell(nST,1);
    sliceFlags=false(nImgs,nST);
    
    for iSlice=1:nST
        sliceFlags(iSlice:nST:nImgs,iSlice)=true;
    end;

refs=GaussFilt(randn(n0,n0,nRefs),.1,nRefs>1);  % random starting refs

%%
iter=0;
%%
while iter<nIters
    %%
    iter=iter+1
    if iter<3
        ri.sigmaC=.2;  % tight prior on the first iterations.
        ri.sigmaG=.2;
    end;
    tic
    %% ------------------------ Do the EM step ---------------------------
    % ------------------------------------------------------------------
    if useParFor
            disp('parallel EM:');
        parfor iSlice=1:nST
            [si1{iSlice},uImg1]=rsStackSplit(sliceFlags(:,iSlice),si,uImgs);
            mo1=moi;
            mo1.imgAmps=moi.imgAmps(si1{iSlice}.origParticle);
            if numel(moi.pRefs)>1
                mo1.pRefs=moi.pRefs(si1{iSlice}.origParticle,:);
            end;
            disp(['EM step ' num2str(iSlice)]);
            [classMeans1{iSlice}, classNorms1{iSlice},...
                ro1(iSlice), imgBestMatch1{iSlice}]...
                =reEMStep22(uImg1,refs,si1{iSlice},ri,mo1);
        end;        
    else
        for iSlice=1:nST
            [si1{iSlice},uImg1]=rsStackSplit(sliceFlags(:,iSlice),si,uImgs);
            mo1=moi;
            mo1.imgAmps=moi.imgAmps(si1{iSlice}.origParticle);
            if numel(moi.pRefs)>1
                mo1.pRefs=moi.pRefs(si1{iSlice}.origParticle,:);
            end;            
            disp(['EM step ' num2str(iSlice)]);
            [classMeans1{iSlice}, classNorms1{iSlice},...
                ro1(iSlice), imgBestMatch1{iSlice}]...
                =reEMStep22(uImg1,refs,si1{iSlice},ri,mo1);

%             
%             
%             
%             [qm1, qn1,qr1, qx1]...
%                 =reEMStep22(uImg1,refs,si1{iSlice},ri,mo1);
%             classMeans1{iSlice}=qm1;
%             classNorms1{iSlice}=qn1;
%             imgBestMatch1{iSlice}=qx1;
%             ro1{iSlice}=qr1;
            
            
        end;
    end;
    toc
    % ------------------------------------------------------------------
    % ------------------------------------------------------------------
            
    %%
    
    imgBestMatch=zeros(n0,n0,nImgs,'single');
    
    roi=ro1(1);
    roi.imgAmps=zeros(1,nImgs,'single');
    roi.pTransUnrot=zeros(nt,nt,nImgs,'single');
    roi.pTrans=zeros(nt,nt,nImgs,'single');
    roi.varNs=zeros(1,nImgs,'single');
    roi.logPX=zeros(1,nImgs,'single');
    roi.pRefs=zeros(nRefs,1,nImgs,'single');
    roi.pAlphas=zeros(nAlphas,nImgs,'single');
    roi.imgClass=zeros(1,nImgs,'single');
    roi.imgTA=zeros(3,nImgs,'single');
    
    for iSlice=1:nST
        inds=si1{iSlice}.origParticle;
        %  Merge the accumulated parameters
        roi.imgAmps(1,inds)=ro1(iSlice).imgAmps;
        roi.pTransUnrot(:,:,inds)=ro1(iSlice).pTransUnrot;
        roi.pTrans(:,:,inds)=ro1(iSlice).pTrans;
        roi.varNs(1,inds)=ro1(iSlice).varNs;
        roi.logPX(1,inds)=ro1(iSlice).logPX;
        roi.pRefs(:,:,inds)=ro1(iSlice).pRefs;
        roi.pAlphas(:,inds)=ro1(iSlice).pAlphas;
        roi.imgClass(1,inds)=ro1(iSlice).imgClass;
        roi.imgTA(:,inds)=ro1(iSlice).imgTA;
        
        imgBestMatch(:,:,inds)=imgBestMatch1{iSlice};
        
    end;
    
    
    % accumulate the fourier projections
    
    
    %%
        for iSlice=1:nSlices;
            if iSlice==1
                classMeans=classMeans1{1};
                classNorms=classNorms1{1};
            else
                classMeans=classMeans+classMeans1{iSlice};
                classNorms=classNorms+classNorms1{iSlice};
            end;
        end;
        classNorms=circshift(classNorms,ceil([n0/2 n0/2 0]));  %2D ifftshift
        
        reShowLatentVars(imgs,refs,ri,roi,iter,[4 10]);
            drawnow;
    %% ----------------Do the Fourier reconstruction----------------
    
    disp('Reconstructions');
    k=1;
    newRefs=real(ifft2(fft2(classMeans)./(classNorms+k)));
    
    % Assign new model values
    moi.sigmaN=sqrt(mean(roi.varNs));
    moi.a=median(roi.imgAmps);
    moi.imgAmps=0*roi.imgAmps+median(roi.imgAmps);
    moi.pRefs=shiftdim(roi.pRefs,2);
    refs=newRefs;
    
    moi
    %
    if flagWriteClasses
            outName=sprintf('%sCls%di.mrc',outPath,iter);
            disp(outName);
    end;
    
    %
    
    figure(10);
    %     subplot(nr,nc,3);
    %     title(['Iteration ' num2str(iter)]);
    if flagSaveFigures
        set(gcf,'paperpositionmode','auto');
        outName=sprintf('%sImgsIter%d.jpg',outPath,iter);
        print('-djpeg','-r200',outName);
        save([outPath 'fsc.mat'], 'fsc');
    end;
    
  if doShowParts  
        figure(1);
        ImagicDisplay2(newRefs,2);
        drawnow;
  end;
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
alphasRot=ri.alphas(mxa)';
[~,mxr]=max(roi.pRefs,[],1);
mxr=mxr(:);
cls=cell(nRefs,1);
for i=1:nRefs
    cls{i}=single(find(mxr==i));
end;

% Make the rotated and aligned stack
% rImgs=zeros(n0,n0,nImgs,'single');
nCls=zeros(nRefs,1);
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
        clsMeans(:,:,j)=clsMeans(:,:,j)/nCls(j);
    end;
end;
if doShowParts
figure(3);
ImagicDisplay2(clsMeans,2);
end;

save([siName(1:end-4) '_classes.mat'],'newRefs', 'clsMeans', 'cls', 'roi', 'ri');

