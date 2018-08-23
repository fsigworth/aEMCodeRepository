% reEMReconstruct.m

% Parameters

doSimulateData=0;
doShowPriors=0;
nIters=8;
flagSaveFigures=false;
flagWriteVolumes=false;
outPath='/Users/fred/EMWork/Simulations/Kv/Reconstructions/';

% Working image sizes
nDs=48;
nCrop=40;  % final working size

symmetry=4;

defaultMbnOffsetA=-70;  % position of membrane relative to RSO particle center
% mapName='3KG2mapsub2.9A.mat';
mapName='KvMap.mat';
% Angle and shift search
angleSteps=[5 5 10];    % step sizes in degrees
% angleSteps=[5 10 20];    % step sizes in degrees
angleLimits=[-20 30 360/symmetry];  % starting alpha and beta; range of gamma.
doIsoFlip=true;
doNormalize=1;
maxShift=3; % Maximum translation searched

fscFlag=1;  % 0: compare with crystal; 1: normal FSC; 2: gold standard
removeShifts=0; % Remove the shifts on starting iteration 3 for efficiency
useEM1=false;  % Use the FFT-based EM step instead of the fast EM2 step
doABReconstruction=0;  % Use Alp's adaptive basis reconstruction

startingRes=40;%%%%%

s0=.003;  % reference amplitude factor

sigmaC=2;
sigmaG=2;

volMaskRadius=12;
volMaskHeight=15;


imgAmp=1;
nVols=1;
sigmaN=1;
useParRotations=1;  % Use parallel processing for the image rotations


% ----------------Get the stack files---------------
if doSimulateData
    
    sm.mDefoci=1;   %Defocus for each micrograph, nMics x 1.
    sm.miIndex=1;   %Micrograph number for each image
    sm.angles=[0 50 0];
    sm.rocks=[5 0];
    sm.bobs=0;
    sm.clicks=[0 0];
    sm.isos=0;
    sm.pixA=4;
    sm.vesR=200/sm.pixA;
    sm.mbnOffset=0;
    sm.sigmaN=1;
    sm.imgAmp=.03;
    sm.n=64;
    sm.mapName='KvMap.mat';
    
    [allImgs,si]=reMakeSimStackSub(sm);
    
    
else
    
    
    [fname, pa]=uigetfile('*si.mat','Select si files','multiselect','on');
    if isnumeric(pa)  % user clicked Cancel
        return
    end;
    if ~iscell(fname)
        fname={fname};
    end;
    cd(pa);
    
    
    
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
    
    si=allSi;
    
end;  % if flagSimulatedData

if ~isfield(si,'mbnOffset')
    si.mbnOffset=defaultMbnOffsetA/si.pixA;
    disp(['Using the default mbn offset: ' num2str(si.mbnOffset*si.pixA) ' A']);
end;


fc=si.pixA/startingRes;  % lowpass for starting volume


n0=size(allImgs,1);
nImgs=size(allImgs,3);

if doNormalize
    % Make an annulus for computing mean and variance
    %%%%% hardwired??
    ringo=45;
    ringi=35;
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

% First, get the first reference volumes.
origVols=arGetRefVolumes(si.pixA,n0,mapName,nVols);
refVols=zeros(size(origVols),'single');
for iVol=1:nVols
    v=s0*SharpFilt(origVols(:,:,:,iVol),fc,fc/n0);
    refVols(:,:,:,iVol)=v;
end;
refVols=repmat(refVols,1,1,1,1,2);  % copy for twin reconstructions


%  We now have imgs and the si struct initialized.
figure(1);
ImagicDisplay1(imgs);
drawnow;

% Set up the model info
moi.sigmaN=sigmaN;
moi.sigmaC=sigmaC;
moi.sigmaG=sigmaG;
moi.b0=0;
moi.pVols=ones(nVols,1)/nVols;
moi.imgAmps=imgAmp*ones(nImgs,1);
moi.a=mean(moi.imgAmps);
% refVols has already been assigned.

fsc=zeros(n0/2-1,nIters,nVols);
for iVol=1:nVols
    figure(4+iVol);
    ShowSections2(refVols(:,:,:,iVol));
    subplot(3,3,9);
    freqs=(0:n0/2-2)/(n0*si.pixA);
    fsc(:,1,iVol)=FSCorr(refVols(:,:,:,iVol),origVols(:,:,:,iVol));
    plot(freqs,fsc(:,:,iVol));
end;

% Set up the run's 2d mask info
ri=struct;
ri.radius=floor(n0/2)-maxShift;
ri.nTrans=2*maxShift+1;
ri.softMask=fuzzymask(n0,2,ri.radius*.9,ri.radius*.2);

%%
ri.volMask=zeros(n0,n0,n0,'single');
ctr=floor(n0/2+1);
for iz=ctr-volMaskHeight:ctr+volMaskHeight
    ri.volMask(:,:,iz)=fuzzymask(n0,2,volMaskRadius,.5);
end;
ri.volMask=GaussFilt(ri.volMask,.05);



% --------Pick angles for the references----------
[ri, refAngles]=reSetRefAngles(angleSteps,angleLimits,doIsoFlip,false,ri);

nRefs=size(refAngles,1);
disp(['nAngles: ' num2str(ri.angleN) '  nRefs: ' num2str(nRefs)]);

tic
uImgs=zeros([n0 n0 nImgs numel(ri.alphas)],'single');

%---------------parfor loop for rotating images----------------
if useParRotations
    disp('Rotating images -parallel');
    parfor i=1:nImgs
        uImgs(:,:,i,:)=rsRotateImage(imgs(:,:,i),-ri.alphas);
    end;
else
    disp('Rotating images');
    for i=1:nImgs
        uImgs(:,:,i,:)=rsRotateImage(imgs(:,:,i),-ri.alphas);
    end;
end;
toc

%% iterate
nTwins=1+single(fscFlag>0);
nRefVols=1+single(fscFlag>1);

outputVols=zeros(n0,n0,n0,nVols,nTwins,nIters);
newVols=zeros(n0,n0,n0,nVols,nTwins);
volSD=zeros(nVols,1);

%%
iter=0;
%%
while iter<nIters
    %%
    iter=iter+1
    %     doABReconstruction= iter>nIters-1;
    doNormalizeVolSD=doABReconstruction;
    
    allProjs=zeros(n0,n0,nRefs,2,nVols,nTwins,'single');
    
    ri.imgAmps=moi.imgAmps;
    ri.n=n0;
    
    if any(iter==removeShifts)  % remove the shifts at this iteration
        disp('Removing shifts');
        if maxShift>0
            [x,y]=ndgrid(-maxShift:maxShift);
            for i=1:nImgs
                img=roi.pTrans(:,:,i);
                ix=x(:)'*img(:);
                iy=y(:)'*img(:);
                p=FourierShift(n0,-[ix iy]);
                if i==1
                    disp([ix iy]);
                end;
                imgs(:,:,i)=real(ifftn(fftn(imgs(:,:,i)).*p));
            end;
        end;
        disp('Rotating images -parallel');
        parfor i=1:nImgs
            uImgs(:,:,i,:)=rsRotateImage(imgs(:,:,i),-ri.alphas);
        end;
        
    end;
    
    %     Pick alternating images
    twinFlags=true(nImgs,1);
    if nTwins>1
        twinFlags(2:nTwins:end,1)=false;
        twinFlags(:,2)=~twinFlags(:,1);
    end;
    
    for i=1:nVols*nTwins
        refVols(:,:,:,i)=refVols(:,:,:,i).*ri.volMask;
    end;
    
    if fscFlag==1  % conventional fsc
        refVols=repmat(mean(refVols,5),1,1,1,1,2); % force the two refVols to be identical.
    end;
    
    %     Normalize the ref amplitudes
    for iVol=1:nVols
        v=mean(refVols(:,:,:,iVol,:),5);  % mean over twins
        vSigma=std(v(:));
        if doNormalizeVolSD && iter>1
            refVols(:,:,:,iVol,:)=refVols(:,:,:,iVol,:)*volSD(iVol)/vSigma;
        end;
        volSD(iVol)=vSigma;
    end;
    
    tic
    for iTwin=1:nTwins
        %         disp('Making refs');
        refs=reMakeTemplates(refVols(:,:,:,:,iTwin),refAngles);  % refs(x,y,iRef,iVol)
        figure(2);
        ImagicDisplay1(refs);
        drawnow;
        
        ri.imgActive=twinFlags(:,iTwin);
        ri.radius=floor(n0/2)-maxShift;
        ri.nTrans=2*maxShift+1;
        ri.pVols=moi.pVols;
        ri.sigmaC=moi.sigmaC;
        ri.sigmaG=moi.sigmaG;
        ri.sigmaN=moi.sigmaN;
        ri.symmetry=symmetry;
        ri.accumClassMeans=false;
        
        %% ------------------------ Do the EM step ---------------------------
        % ------------------------------------------------------------------
        disp(['EM step ' num2str(iTwin)]);
        if useEM1
            [classMeans, classNorms, ro1, pTrans1, pAlphas1, pRefs1]=reEMStep1(imgs,refs,si,ri);
        else
            if doSimulateData  % drop out after looking at the first image.
                [classMeans, classNorms, ro1, rawMeans, rawNorms, imgBestMatch1,rawLPs1]...
                    =reEMStep22(uImgs,refs,si,ri);
                rawLProbs=rawLPs1;
                figure(3); SetGrayscale;
                reShowDebug;
                return
            else
                [classMeans, classNorms, ro1, rawMeans, rawNorms, imgBestMatch1]...
                    =reEMStep22(uImgs,refs,si,ri);
            end;
        end;
        
        % ------------------------------------------------------------------
        % ------------------------------------------------------------------
        %  Combine the accumulated parameters
        if iTwin==1
            roi=ro1;
            imgBestMatch=imgBestMatch1;
            %             rawLProbs=rawLPs1;
        else
            roi=AddStructFields(roi,ro1);
            imgBestMatch=imgBestMatch+imgBestMatch1;
            %             rawLProbs=AddStructFields(rawLProbs,rawLPs1);
        end;
        
        
        if doABReconstruction
            nMi=size(si.ctfs,3);
            sumMeans=sum(rawMeans,5);  % sum over all the ctfs
            sumCtfs=zeros(n0,n0,nRefs,nVols,'single');
            %             repmat(shiftdim(classNorms,-2),n0,n0,1,1);
            normNorms=sum(rawNorms,3)+1;  % norm=1 for 1 image contributing.
            for j=1:nRefs
                for k=1:nVols
                    for i=1:nMi
                        sumCtfs(:,:,j,k)=sumCtfs(:,:,j,k)+si.ctfs(:,:,i)*rawNorms(j,k,i);
                    end;
                end;
            end;
            %             We normalize ctfs to be the same as one ctf.  We normalize
            %             class means to be ~sqrt(nImgsInClass)
            for j=1:nRefs;
                for k=1:nVols
                    sumCtfs(:,:,j,k)=sumCtfs(:,:,j,k)/normNorms(j,k);
                    sumMeans(:,:,j,k)=sumMeans(:,:,j,k)/sqrt(normNorms(j,k));
                end;
            end;
            
            disp('abReconstruct');
            for iVol=1:nVols
                %%
                figure(3);
                ImagicDisplay2(sumMeans(:,:,:,iVol));
                drawnow;
                
                newVols(:,:,:,iVol,iTwin)=abReconstruct(...
                    sumMeans(:,:,:,iVol),sumCtfs(:,:,:,iVol),refAngles);
                figure(5);
                ShowSections2(newVols(:,:,:,iVol,iTwin));
                title(['ABReconstruct ' num2str([iVol iTwin])]);
            end;
        else % accumulate the fourier projections
            
            %%
            for iVol=1:nVols
                realNorms=zeros(n0,n0,nRefs,'single');
                nCls=size(classNorms,3);
                for i=1:nCls
                    realNorms(:,:,i)=fftshift(real(ifftn(ifftshift(classNorms(:,:,i,iVol)))));
                end;
                classMeanSym=reSymmetrizeClasses(classMeans,ri);
                figure(3);
                ImagicDisplay2(classMeanSym(:,:,:,iVol));
                drawnow;
                
                allProjs(:,:,:,1,iVol,iTwin)=classMeanSym(:,:,:,iVol);
                %              allProjs(:,:,:,1,iVol,iTwin)=classMeans(:,:,:,iVol);
                allProjs(:,:,:,2,iVol,iTwin)=realNorms;
            end;
        end; % if abReconstruct
    end; % iTwin
    toc
    %
    reShowLatentVars(imgs,refs,ri,roi,iter,[4 10]);
    drawnow;
    
    %    tic
    
    if ~doABReconstruction
        %% ----------------Do the Fourier reconstruction----------------
        
        disp('Reconstructions');
        allFVols=reFourierInsertion(allProjs,refAngles,ri.symmetry);
        %
        %     Normalization of volumes
        k=100;
        for iTwin=1:nTwins
            for iVol=1:nVols
                newVols(:,:,:,iVol,iTwin)=rsNormalizeReconstruction(allFVols(1,iVol,iTwin),allFVols(2,iVol,iTwin),k);
            end;
        end;
        toc
        
        figure(6);  % Show the 3D Fourier norms
        ShowSections2(allFVols(2).PadFT);
        subplot(3,3,1);
        title('Fourier norms');
    end;
    
    % Display the volumes and FSCs
    freqs=(0:n0/2-2)'/(n0*si.pixA);  % frequencies for fscs
    volumeProbs=mean(roi.pVols,2)
    % Display volumes and fscs
    for iVol=1:nVols
        
        
        figure(4+iVol);
        ShowSections2(mean(newVols(:,:,:,iVol,:),5));
        
        subplot(3,3,1);
        title(['p(vol)= ' num2str(volumeProbs(iVol))]);
        
        subplot(3,3,9);
        if fscFlag>0
            fsc(:,iter+1,iVol)=FSCorr(newVols(:,:,:,iVol,1),newVols(:,:,:,iVol,2));
            target=0.143;
        else
            fsc(:,iter+1,iVol)=FSCorr(newVols(:,:,:,iVol),origVols(:,:,:,iVol));
            target=0.5;
        end;
        plot(freqs,[fsc(:,:,iVol) 0*freqs+target 0*freqs]);
        hold on;
        plot(freqs,fsc(:,iter+1,iVol),'k.-','markersize',8);
        hold off;
        
        title('FSC');
        outputVols(:,:,:,iVol,iTwin,iter)=newVols(:,:,:,iVol,iTwin);
    end;
    
    % Assign new model values
    refVols=newVols;
    moi.sigmaN=sqrt(mean(roi.varNs));
    moi.a=median(roi.imgAmps);
    %     moi.sigmaC=sqrt(mean(roi.sigmaC));
    % sigmaG remains
    %     moi.b0=mean(roi.ePars(:,4));
    %     moi.pVols=mean(roi.pVols);
    moi.imgAmps=0*roi.imgAmps+median(roi.imgAmps);
    
    moi
    %
    if flagWriteVolumes
        for iVol=1:nVols
            outName=sprintf('%svol%di%d.mrc',outPath,iVol,iter);
            disp(outName);
            WriteMRC(newVols(:,:,:,iVol),si.pixA,outName);
        end;
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
    
    
    figure(1);
    ImagicDisplay2(imgBestMatch);
    %     q=reshape(roi.pRefs,nRefs*nVols,nImgs);
    % [vals, roi.classInds]=max(q);
    % roi.classInds=roi.classInds(:);
end; % while iter
