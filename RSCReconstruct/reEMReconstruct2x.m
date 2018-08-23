% reEMReconstruct1.m

% Parameters
fl=struct;
fl.makeFakeData=false;
fl.doShowPriors=false;
fl.useAltImgs=true;  % use the unsubtracted stack for alternative reconstruction.
fl.removeRings=true;

fl.saveFigures=true;
fl.writeVolumes=true;
fl.useParFor=true;

ri=struct;
ri.flag=fl;
ri.nJobs=1;
ri.nSlices=4;
ri.iTwin=0;  % undefined; values are 1 and 2

pr.transThreshold=1e-5;

% ri.fscFlag=MyInput('fscFlag',-1);  % 0: compare with crystal; 1: normal FSC; 2: gold standard
%             -1: plain reconstruct from odd images; -2: even
ri.startIter=1;
ri.nIters=10;

ri.stackNames={'sq10_350p128tsi.mat'};
ri.stackPath='Stack/';

ri.outPath='Reconstructions/Recon40a/'; % relative to the root path
ri.outName='141106';
ri.kFactor=1;


% Working image sizes
% ri.nDs=48;
ri.nDs=128;
% nCrop=40;  % final working size
ri.nCrop=96;

volMaskRadius=.25*ri.nDs;
volMaskHeight=round(.35*ri.nDs);
ri.volMask=zeros(n0,n0,n0,'single');
ctr=floor(n0/2+1);
for iz=ctr-volMaskHeight:ctr+volMaskHeight
    ri.volMask(:,:,iz)=fuzzymask(n0,2,volMaskRadius,.5);
end;
ri.volMask=GaussFilt(ri.volMask,.05);

fscMaskRadius=.33*ri.nDs;
fscMaskWidth=.08*ri.nDs;
ri.fscMask=fuzzymask(n0,3,fscMaskRadius,fscMaskWidth);

ri.symmetry=4;

% % defaultMbnOffsetA=-70;  % position of membrane relative to RSO particle center
ri.startMapName='KvMap.mat';
ri.startRes=60;

% Angle and shift search
angleSteps=[5 5 10]
% angleSteps=[5 5 10];    % step sizes in degrees
% angleSteps=[3 3 5]    % step sizes in degrees
angleLimits=[10 30 360/symmetry];  % starting alpha and beta; range of gamma.
% isos=[0 1];
isos=1;
ri=reSetRefAngles(angleSteps,angleLimits,isos,false,ri);

doNormalize=1;
maxShift=2; % Maximum translation searched
% maxShift=3; % Maximum translation searched
nt=2*maxShift+1;

removeShifts=0; % Remove the shifts on this siteration for efficiency

startingRes=50; %%%%%%%%%%%%

s0=.1;  % reference amplitude factor
sim.imgScale=.015;
sim.imgScale=.05;
sim.pIso=0;
sim.n=64;
sim.nMicrographs=5;
sim.sigmaN=1;
sim.sigmaC=2;
sim.sigmaA=5;
sim.angleSteps=angleSteps;
sim.angleLimits=angleLimits;
sim.nMi=1;

sigmaC=2;
sigmaG=2;

imgAmp=1;
nVols=1;
sigmaN=10;  % start with a high value to favor priors
% sigmaN=1;  % start with a high value to favor
simSigmaN=1;


% -----See if we're restarting-----
iPhase=max(-ri.fscFlag,0);
outNameAll=sprintf('%s%sOutputVols%dt.mat',ri.outPath,ri.prefix,iPhase);
doRestart=exist(outNameAll,'file');
if doRestart
    disp(['Found the restart file:' outNameAll]);
    load(outNameAll);  % get the outputVols
end;

%% make fake data?
if makeFakeData
    disp('Making simulated data');
    [si0,imgs0]=reMakeFakeData(sim,mapName);
else  % if names and stackPath already exist, use them.
    if numel(ri.names)<1
        ri.names=[];
    end;

    % ----------------Get the stack files---------------
    %     ri.names={'sq10_180_Jun25_20.18.50Ctsi.mat'};
    %     [si0, imgs0, ri.names]=reLoadStackFiles(struct,ri.names);
    if useAltImgs
        [si0,imgs0,ri.names,ri.stackPath,altImgs0]=reLoadStackFiles(struct,ri.names,ri.stackPath);
    else
        [si0,imgs0,ri.names,ri.stackPath]=reLoadStackFiles(struct,ri.names,ri.stackPath);
    end;
    disp(ri.names);
end;

if ri.fscFlag<0  % pick up odd or even images
    [nm0, nac0]=size(si0.activeFlags);
    si0.activeFlags(:,nac0+1)=si0.activeFlags(:,nac0);
    startInd=max(1,-ri.fscFlag);  %fscFlag=-1 means odd, -2 means even are kept
    si0.activeFlags(startInd:2:nm0,end)=false;
end;
[si,imgs]=rsStackSplit(si0.activeFlags(:,end),si0,imgs0);
[si, imgs]=rsStackDownsample(si,imgs,ri.nDs);
[si, imgs]=rsStackCrop(si,imgs,nCrop);

%
if useAltImgs
    [six,altImgs0]=rsStackSplit(si0.activeFlags(:,end),si0,altImgs0);
    [six, altImgs]=rsStackDownsample(six,altImgs0,ri.nDs);
    [six, altImgs]=rsStackCrop(six,altImgs,nCrop);
else
    altImgs=[];
end;
clear si0 % save space.

if flagShowImgs
    disp('Showing images');
    figure(1);
    ImagicDisplay3(imgs);
    drawnow;
end;

if ~exist(ri.outPath,'dir')
    mkdir(ri.outPath);
end;

n0=size(imgs,1)
nImgs=size(imgs,3)

rotMask=fuzzymask(n0,2,n0/2,n0/20);

if ~exist(ri.outPath,'dir')
    disp(['Making the directory ' ri.outPath]);
    mkdir(ri.outPath);
end;


%% --------Initialize the model ---------

% Get the first reference volumes.
[origVols,mbnOffsetA,volSD]=arGetRefVolumes(si.pixA,n0,mapName,nVols);
refVols=zeros(size(origVols),'single');
fc=si.pixA/startingRes;  % lowpass for starting volume
volSD=zeros(nVols,1);
for iVol=1:nVols
    v=SharpFilt(origVols(:,:,:,iVol),fc,fc/n0);
    refVols(:,:,:,iVol)=v;
    volSD(iVol)=sqrt(v(:)'*v(:)/numel(v))
end;

nTwins=1+single(ri.fscFlag>0);
refVols=repmat(refVols,1,1,1,1,nTwins);  % copy for twin reconstructions

% Set up the other model pars
moi.sigmaN=sigmaN;
moi.sigmaC=sigmaC;
moi.sigmaG=sigmaG;
moi.b0=0;
moi.pVols=ones(nVols,1)/nVols;
moi.imgAmps=imgAmp*ones(nImgs,1);
moi.a=mean(moi.imgAmps);
moi.pRefs=1;
moi.activeTrans=true(nt,nt,nImgs);

fsc=zeros(n0/2-1,ri.nIters,nVols);
for iVol=1:nVols
    figure(4+iVol);
    ShowSections2(refVols(:,:,:,iVol));
    subplot(3,3,9);
    freqs=(0:n0/2-2)/(n0*si.pixA);
    fsc(:,1,iVol)=FSCorr(refVols(:,:,:,iVol),origVols(:,:,:,iVol));
    plot(freqs,fsc(:,:,iVol));
end;

% ----------Set up the Run Info structure-----------
ri.radius=floor(n0/2)-maxShift;
ri.nTrans=nt;
ri.softMask=fuzzymask(n0,2,ri.radius*.9,ri.radius*.2);
ri.symmetry=symmetry;
ri.useAltImgs=0;  % start off using the normal images.


% Masks

% Pick angles for the references
[ri, refAngles]=reSetRefAngles(angleSteps,angleLimits,isos,false,ri);

nRefs=size(refAngles,1);
disp(['nAngles: ' num2str(ri.angleN) '  nRefs: ' num2str(nRefs)]);
disp(['shift range: ' sprintf('%d x %d',nt,nt)]);


%--------------Set up slicing of the dataset---------------
nAlphas=numel(ri.alphas);
nAI=numel(ri.alphasI);
newVols=zeros(n0,n0,n0,nVols,nTwins);

nST=nSlices*nTwins;
nSliceImgs=ceil(nImgs/nST);

ro1=cell(nST,1);
imgBestMatch1=cell(nST,1);
classMeans1=cell(nST,1);
classNorms1=cell(nST,1);
si1=cell(nST,1);

%     Slices are identified by booleans.  Alternating images are assigned
%     to iTwin=1,2; blocks of these are then assigned to slices.
sliceFlags=false(nImgs,nTwins,nSlices);
numInSlice=ceil(nImgs/nST)*nTwins;
for is=1:nSlices
    sStart=(is-1)*numInSlice;
    sEnd=min(is*numInSlice,nImgs);
    for it=1:nTwins
        sliceFlags(sStart+it:nTwins:sEnd,it,is)=true;
    end;
end;

ring.radii=[-30 -20 -10 0 10 20 30]/si.pixA;
ring.widths=ones(size(ring.radii))*max(1,14/si.pixA);
ring.fitMask=(0.8*ri.softMask+0.2);
ringFits=cell(nST,1);
for i=1:nST
    ringFits{i}=single(0);
end;
figure(10);
%     Make the figure 800 x 800 pixels
pos=get(gcf,'position');
pos(3:4)=[840 750];
set(gcf,'position',pos);

% Are we restarting?
if doRestart
    iter=numel(outputVols);
    fprintf(' restarting after iteration %d\n',iter);
    newVols=outputVols{iter};
    refVols=newVols;
    roi=rois{iter};
    moi=mois{iter,1};
    iter=iter+1;
else
    iter=1;
    mois={moi};
    rois=cell(0,1);
    outputVols=cell(0,1);
end;
%%
while iter<=ri.nIters
    iteration=iter
    ri.useAltImgs=(iter==ri.nIters) && useAltImgs;  % turn on use of alt images
    % at the last iteration
    %
    allProjs=zeros(n0,n0,nRefs,2,nVols,nTwins,'single'); % 2 is sums/norms
    
    %    multiply by the volume mask
    refVols=refVols.*repmat(ri.volMask,1,1,1,nVols,nTwins);
    %     Set up twin volumes if necessary
    if ri.fscFlag==1  % conventional fsc
        refVols=repmat(mean(refVols,5),1,1,1,1,2); % force the two refVols to be identical.
    end;
    
    %     Normalize the ref amplitudes
    refVolsNorm=refVols;
    for iVol=1:nVols
        v=mean(refVols(:,:,:,iVol,:),5);  % mean over twins
        vSigma=sqrt(v(:)'*v(:)/numel(v))
        %         vSigma=volSD(iVol); %%%%%%
        refVolsNorm(:,:,:,iVol,:)=refVols(:,:,:,iVol,:)*volSD(iVol)/vSigma;
    end;
    figure(4);
    ShowSections2(refVolsNorm(:,:,:,1));
    title(['Normalized refs ' num2str(iter)]);
    drawnow;
    
    tic
    
    disp('Making refs');
    refs=s0*reMakeTemplates(refVolsNorm,refAngles);  % refs(x,y,iRef,iVol,iTwin)
    
    %% ------------------------ Do the EM step ---------------------------
    % ------------------------------------------------------------------
    tic
    if useParFor
        disp('parallel EM');
        parfor( is=1:nST,nST)
            [si1{is},img1,altImg1]=rsStackSplit(sliceFlags(:,is),si,imgs,altImgs);
            iTw=mod(is,nTwins)+1;  % twin index
            mo1=reModelSplit(moi,si1{is}.pastParticle);
            [classMeans1{is}, classNorms1{is},ro1{is}, imgBestMatch1{is}]...
                =reEMStep23(img1-ringFits{is},refs(:,:,:,:,iTw),si1{is},ri,mo1,altImg1);
            if removeRings
                mRefs=reMakeMatchedRefs(refs,si1{is},ri,ro1{is});
                ringFits{is}=reFitVesicleRings(si1{is},img1-mRefs,ring.fitMask,ring.radii,ring.widths);
            end;
        end;
    else
        disp('EM');
        for is=1:nST
            disp(is);
            %             [si1{is},uImg1,uAltImg1]=rsStackSplit(sliceFlags(:,is),si,uImgs,uAltImgs);
            %             iTw=(is>nSlices)+1;  % twin index
            %             mo1=reModelSplit(moi,si1{is}.pastParticle);
            %             [classMeans1{is}, classNorms1{is},ro1{is}, imgBestMatch1{is}]...
            %                 =reEMStep22(uImg1,refs(:,:,:,:,iTw),si1{is},ri,mo1,uAltImg1);
            [si1{is},img1,altImg1]=rsStackSplit(sliceFlags(:,is),si,imgs,altImgs);
            iTw=mod(is,nTwins)+1;  % twin index
            mo1=reModelSplit(moi,si1{is}.pastParticle);
            [classMeans1{is}, classNorms1{is},ro1{is}, imgBestMatch1{is}]...
                =reEMStep23(img1-ringFits{is},refs(:,:,:,:,iTw),si1{is},ri,mo1,altImg1);
            if removeRings
                mRefs=reMakeMatchedRefs(refs,si1{is},ri,ro1{is});
                ringFits{is}=reFitVesicleRings(si1{is},img1-mRefs,ring.fitMask,ring.radii,ring.widths);
            end;
        end;
    end;
    toc
    % ------------------------------------------------------------------
    % ------------------------------------------------------------------
    
    %% -----Accumulate slices----------
    
    roi=struct;
    roi.imgAmps=zeros(1,nImgs,'single');
    %     roi.pTransUnrot=zeros(nt,nt,nImgs,'single');
    roi.pTrans=zeros(nt,nt,nImgs,'single');
    roi.varNs=zeros(1,nImgs,'single');
    roi.pVols=zeros(nVols,nImgs,'single');
    roi.logPX=zeros(1,nImgs,'single');
    roi.pRefs=zeros(nRefs,nVols,nImgs,'single');
    roi.pAlphas=zeros(nAI,nImgs,'single');
    roi.imgClass=zeros(1,nImgs,'single');
    roi.imgTA=zeros(3,nImgs,'single');
    imgBestMatch=zeros(n0,n0,nImgs,'single');
    
    for is=1:nST
        inds=si1{is}.pastParticle;
        roi.imgAmps(1,inds)=ro1{is}.imgAmps;
        %         roi.pTransUnrot(:,:,inds)=ro1{is}.pTransUnrot;
        roi.pTrans(:,:,inds)=ro1{is}.pTrans;
        roi.varNs(1,inds)=ro1{is}.varNs;
        roi.pVols(:,inds)=ro1{is}.pVols;
        roi.logPX(1,inds)=ro1{is}.logPX;
        roi.pRefs(:,:,inds)=ro1{is}.pRefs;
        roi.pAlphas(:,inds)=ro1{is}.pAlphas;
        roi.imgClass(1,inds)=ro1{is}.imgClass;
        roi.imgTA(:,inds)=ro1{is}.imgTA;
        imgBestMatch(:,:,inds)=imgBestMatch1{is};
    end;
    
    
    % ---------Accumulate class means-----------
    for iTwin=1:nTwins
        iSOffset=nSlices*(iTwin-1);
        for is=iSOffset+1:iSOffset+nSlices;
            if is==iSOffset+1
                classMeans=classMeans1{1};
                classNorms=classNorms1{1};
            else
                classMeans=classMeans+classMeans1{is};
                classNorms=classNorms+classNorms1{is};
            end;
        end;
        for iVol=1:nVols
            realNorms=zeros(n0,n0,nRefs,'single');
            nCls=size(classNorms,3);
            for i=1:nCls
                realNorms(:,:,i)=fftshift(real(ifftn(ifftshift(classNorms(:,:,i,iVol)))));
            end;
            %             classMeanSym=reSymmetrizeClasses(classMeans,ri);
            allProjs(:,:,:,1,iVol,iTwin)=classMeans(:,:,:,iVol);
            allProjs(:,:,:,2,iVol,iTwin)=realNorms;
        end;
    end; % iTwin
    if flagShowImgs
        figure(2);
        ImagicDisplay3(squeeze(sum(allProjs(:,:,:,1,:,:),6)));
    end;
    %
    
    reShowLatentVars(imgs,refs,ri,roi,iter,[7 10],['iter ' num2str(iter)]);
    drawnow;
    
    
    
    %% ----------------Do the Fourier reconstruction----------------
    disp('Reconstructions');
    tic
    usePar=useParFor && (nVols*nTwins>=4);  % don't use parFor if few volumes
    allFVols=reFourierInsertion(allProjs,refAngles,ri.symmetry,usePar);
    %
    %     Normalization of volumes
%     k=median(roi.varNs)/4;
    k=mean(roi.imgAmps.^2)*ri.kFactor  % mult by fsc/(1-fsc) says Sjors
    for iTwin=1:nTwins
        for iVol=1:nVols
            newVols(:,:,:,iVol,iTwin)=rsNormalizeReconstruction(allFVols(1,iVol,iTwin),allFVols(2,iVol,iTwin),k);
        end;
    end;
    toc
    
    figure(6);  % Show the 3D Fourier norms
    ShowSections2(Crop(allFVols(2).PadFT,n0));
    subplot(3,3,1);
    title('Fourier norms');
    
    % Display the volumes and FSCs
    freqs=(0:n0/2-2)'/(n0*si.pixA);  % frequencies for fscs
    volumeProbs=mean(roi.pVols,2)
%%    % Display volumes and fscs

    for iVol=1:nVols
        figure(4+iVol);
        ShowSections2(fscMask.*mean(newVols(:,:,:,iVol,:),5));
        
        subplot(3,3,1);
        title(['Reconstruction ' num2str(iter)]);
        
        subplot(3,3,9);
        if ri.fscFlag>0
            fsc(:,iter+1,iVol)=FSCorr(fscMask.*newVols(:,:,:,iVol,1),fscMask.*newVols(:,:,:,iVol,2));
            target=0.143;
        else
            fsc(:,iter+1,iVol)=FSCorr(fscMask.*newVols(:,:,:,iVol),fscMask.*origVols(:,:,:,iVol));
            target=0.5;
        end;
        plot(freqs,[fsc(:,:,iVol) 0*freqs+target 0*freqs]);
        hold on;
        plot(freqs,fsc(:,iter+1,iVol),'k.-','markersize',8);
        hold off;
        title('FSC');
        drawnow;
    end;
    %%
    % Assign new model values
    %     updateRefVols=(iter>0);
    %     if updateRefVols
    refVols=newVols;
    %     end;
    moi.sigmaN=sqrt(mean(roi.varNs));
    %     moi.sigmaC=...
    moi.a=median(roi.imgAmps);
%     moi.imgAmps=0*roi.imgAmps+median(roi.imgAmps);
    moi.imgAmps=reSetImageAmplitudes(si,roi.imgAmps);
    moi.activeTrans=reAssignActiveTrans(roi.pTrans,transThreshold);
    %     moi
    fracActiveTrans=sum(moi.activeTrans(:))/numel(moi.activeTrans)
%     fracActiveAlphas=sum(roi.pAlphas(:)>transThreshold)/numel(roi.pAlphas)
%     fracActiveRefs=sum(roi.pRefs(:)>transThreshold)/numel(roi.pRefs)
    mois{iter,1}=moi;
    %

    figure(10);
    
    %     subplot(nr,nc,3);
    %     title(['Iteration ' num2str(iter)]);
    if flagSaveFigures
        iTwin=max(1,-ri.fscFlag);
        set(gcf,'paperpositionmode','auto');
        outName=sprintf('%s%sIter%02dt%d.jpg',ri.outPath,ri.prefix,iter,iTwin);
        print('-djpeg','-r200',outName);
        save([ri.outPath 'fsc.mat'], 'fsc');
    end;
    
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
    else
        corrImgs=imgs;
    end;
      
    if flagWriteVolumes
        outNameAll=sprintf('%s%sOutputVols%dt.mat',ri.outPath,ri.prefix,iPhase);
        for iVol=0:nVols*nTwins
            outName=sprintf('%s%sv%dt%di%02d.mrc',ri.outPath,ri.prefix,iVol,iPhase,iter);
            disp(outName);
            if iVol==0
                if nTwins>1
                    WriteMRC(mean(newVols(:,:,:,1,:),6),si.pixA,outName);
                end;
            else
                WriteMRC(newVols(:,:,:,iVol),si.pixA,outName);
            end;
        end;
        
        outputVols{iter,1}=newVols;
        rois{iter,1}=roi;
        if iter>1  % for previous iterations, save only the mean of pRefs
            rois{iter-1,1}.pRefs=sum(rois{iter-1,1}.pRefs,3);
        end;
        save(outNameAll,'outputVols','mois','rois','si','ri','fsc');
        outNameRings=sprintf('%s%sOutputVols%dt.mat',ri.outPath,ri.prefix,iPhase);
    end;
    

    
    
    
    iter=iter+1;
end; % while iter
