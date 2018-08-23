% Get env variables $EM_WORK_DIR, $EM_NODE_NUM
% In WORK_DIR expect to find ari.mat, vols.mat volsI01M01.mrc etc.,
%  imgs altImgs, si01
%  mois.mat rings.mrc and /scratch
% scratch contains si01, si02...
% ari has full path to tsi.mat, tsi.mrc, tusi.mrc
% ari.numNodes
% Write scratch/availN01.mat (node 01 is alive)
% Write scratch/roN01.mat, clsN01.mrc,

% Find out if we're a node for a parallel job
numNodes=str2double(getenv('PBS_NUM_NODES'));
if isnan(numNodes)
    numNodes=1;
    iNode=1;
    workDir='';
    scratchDir='scratch/';
    jobName='EMReconNode';
else
    iNode=str2double(getenv('PBS_NODE_INDEX'));
    workDir=getenv('EM_WORK_DIR');
    jobName=getenv('EM_JOBNAME');
    scratchDir=getenv('EM_SCRATCH_DIR');
end;

if iNode<2 % we are running ourselves
    % reEMReconstruct2.m
    
    % Parameters
    nSets=numNodes;
    makeFakeData=0;
    doShowPriors=0;
    useAltImgs=0;  % use the unsubtracted stack for alternative reconstruction.
    removeRings=1;
    nIters=10;
    flagSaveFigures=1;
    flagWriteVolumes=1;
    flagSkipFileSelection=1;
    flagShowImgs=0;
    % outPath='/Users/fred/EMWork/Simulations/Kv/Reconstructions/';
    outPath='';
    useParFor=1;
    nSlices=4;
    prefix='KvRings1';
    if ~exist('names','var')
        names=cell(0,0);
    end;
    if ~exist('siPath','var')
        siPath='';
    end;
    % names='sq10_2_350n64tsi.mat';
    
    % Working image sizes
    nDs=48;
    nCrop=40;  % final working size
    volMaskRadius=.25*nDs;
    volMaskHeight=round(.35*nDs);
    
    symmetry=4;
    
    defaultMbnOffsetA=-70;  % position of membrane relative to RSO particle center
    mapName='KvMap.mat';
    
    % Angle and shift search
    angleSteps=[5 5 10];    % step sizes in degrees
    % angleSteps=[2 5 10];    % step sizes in degrees
    angleLimits=[-10 20 360/symmetry];  % starting alpha and beta; range of gamma.
    isos=[0 1];
    doNormalize=1;
    maxShift=4; % Maximum translation searched
    maxShift=3; % Maximum translation searched
    nt=2*maxShift+1;
    
    fscType=1;  % 0: compare with crystal; 1: normal FSC; 2: gold standard
    removeShifts=0; % Remove the shifts on this siteration for efficiency
    
    startingRes=10; %%%%%%%%%%%%
    
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
    
    %% make fake data?
    if makeFakeData
        disp('Making simulated data');
        [si0,imgs0]=reMakeFakeData(sim,mapName);
    else
        % ----------------Get the stack files---------------
        %     names={'sq10_180_Jun25_20.18.50Ctsi.mat'};
        %     [si0, imgs0, names]=reLoadStackFiles(struct,names);
        if useAltImgs
            [si0,imgs0,names,siPath,altImgs0]=reLoadStackFiles(struct,names,siPath);
        else
            [si0,imgs0,names,siPath]=reLoadStackFiles(struct,names,siPath);
        end;
        disp(names);
    end;
    
    [si,imgs]=rsStackSplit(si0.activeFlags(:,end),si0,imgs0);
    [si,imgs]=rsStackDownsample(si,imgs,nDs);
    [si,imgs]=rsStackCrop(si,imgs,nCrop);
    WriteMRC(imgs,si.pixA,[jobName '_imgs.mrc']);
    %
    if useAltImgs
        [six,altImgs0]=rsStackSplit(si0.activeFlags(:,end),si0,altImgs0);
        [six, altImgs]=rsStackDownsample(six,altImgs0,nDs);
        [six, altImgs]=rsStackCrop(six,altImgs,nCrop);
    else
        altImgs=[];
    end;
    
    if flagShowImgs
        disp('Showing images');
        figure(1);
        ImagicDisplay3(imgs);
        drawnow;
    end;
    
    n0=size(imgs,1)
    nImgs=size(imgs,3)
    
    rotMask=fuzzymask(n0,2,n0/2,n0/20);
    
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
        for iTwin=1:2
            volName=[prefix '_vol' sprintf('_i%02v%02t%02.mrc',iVol,0,iTwin)];
            WriteMRC(v,si.pixA,volName);
        end;
    end;
    refVols=repmat(refVols,1,1,1,1,2);  % copy for twin reconstructions
    
    % Set up the other model pars
    moi.sigmaN=sigmaN;
    moi.sigmaC=sigmaC;
    moi.sigmaG=sigmaG;
    moi.b0=0;
    moi.pVols=ones(nVols,1)/nVols;
    moi.imgAmps=imgAmp*ones(nImgs,1);
    moi.a=mean(moi.imgAmps);
    moi.pRefs=1;
    
    fsc=zeros(n0/2-1,nIters,nVols);
    for iVol=1:nVols
        figure(4+iVol);
        ShowSections2(refVols(:,:,:,iVol));
        subplot(3,3,9);
        freqs=(0:n0/2-2)/(n0*si.pixA);
        fsc(:,1,iVol)=FSCorr(refVols(:,:,:,iVol),origVols(:,:,:,iVol));
        plot(freqs,fsc(:,:,iVol));
    end;
    
    % ----------Set up the Run Info structure-----------
    
    
    ri=struct;
    ri.radius=floor(n0/2)-maxShift;
    ri.nTrans=nt;
    ri.softMask=fuzzymask(n0,2,ri.radius*.9,ri.radius*.2);
    ri.symmetry=symmetry;
    ri.useAltImgs=0;  % start off using the normal images.
    ri.fscFlag=fscType;
    ri.scratchDir='scratch/';
    ri.prefix=jobName;
    if ~exist(ri.scratchDir,'dir')
        mkdir(ri.scratchDir);
    end;
    
    ring=struct;
    ring.radii=[-30 -20 -10 0 10 20 30]/si.pixA;
    ring.widths=ones(size(ring.radii))*max(1,14/si.pixA);
    ring.fitMask=(0.8*ri.softMask+0.2);
    ri.ring=ring;

    % Masks
    ri.volMask=zeros(n0,n0,n0,'single');
    ctr=floor(n0/2+1);
    for iz=ctr-volMaskHeight:ctr+volMaskHeight
        ri.volMask(:,:,iz)=fuzzymask(n0,2,volMaskRadius,.5);
    end;
    ri.volMask=GaussFilt(ri.volMask,.05);
    
    % Pick angles for the references
    ri=reSetRefAngles(angleSteps,angleLimits,isos,false,ri);
    
    %--------------Write out files for workers---------------
    numInSet=2*ceil(nImgs/(2*nSets)); % want it to be even
    allSetFlags=false(nImgs,nSets);
    setNames=cell(nSets,1);
    for is=1:nSets
        sStart=(is-1)*numInSet;
        sEnd=min(is*numInSet,nImgs);
        allSetFlags(sStart+1:sEnd,is)=true;
        setNames{is}=sprintf('%02d.mat',is);
        ris=ri;
        ris.setFlags=allSetFlags(:,is);
        save([ri.scratchDir ri.prefix '_ris' setNames{is}],'ris');
        sis=rsStackSplit(ris.setFlags,si);
        save([ri.scratchDir ri.prefix '_sis' setNames{is}],'sis');
        save([ri.scratchDir ri.prefix '_sts' setNames{is}],
    end;
    ri.sliceFlags=allSetFlags;
    save([ri.prefix '_ri.mat'],'ri');
    
    nRefVols=1+single(ri.fscFlag>1);
    nAlphas=numel(ri.alphas);
    nAI=numel(ri.alphasI);
    
    nTwins=1+single(ri.fscFlag>0);
    outputVols=zeros(n0,n0,n0,nVols,nTwins,nIters);
    newVols=zeros(n0,n0,n0,nVols,nTwins);
    
end;  % special code for iNode<2
% At this point we have the following files for node XX
%   EMRecon_risXX.mat
%   EMRecon_sisXX.mat


% -------------Everyone does the following----------------
riName=sprintf('%s%s_ris%02d.mat',scratchDir,jobName,iNode);
siName=sprintf('%s%s_sis%02d.mat',scratchDir,jobName,iNode);
imgsName=[scratchDir,jobName '_imgs.mrc'];

load(riName);  % get ri
load(siName);  % get si

refAngles=reGetAngleList(ri);
    nRefs=size(refAngles,1);
    disp(['nAngles: ' num2str(ri.angleN) '  nRefs: ' num2str(nRefs)]);
    disp(['shift range: ' sprintf('%d x %d',nt,nt)]);

% Load the input files
iSet=max(iNode,1);
sliceName=sprintf('_%02d.mat',iSet);
load([scratchDir jobName '_ris' sliceName]);  % load ris
load([scratchDir jobName '_sis' sliceName]);  % load sis
imgs=ReadMRC([jobName '_imgs.mrc']);
imgs=imgs(:,:,ris.sliceFlags);
sz=size(imgs);
n0=sz(1);
nImgs=sz(3);
if ris.iter>1
    sliceIter=sprintf('_%02di%02d',iSet,ris.iter);
    rings=ReadMRC([ris.prefix sliceIter 'rings.mrc']);
else
    rings=zeros(1,1,nImgs,'single');
end;
refVols=zeros([n0 n0 n0 ri.nVols 2],'single')
for iVol=1:ri.nVols
    for iTwin=1:2
        volName=[prefix '_vol' sprintf('_i%02v%02t%02.mrc',iVol,0,iTwin)];
        refVols(:,:,:,iVol,iTwin)=ReadMRC(volName);
    end;
end;

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

ringFits=cell(nST,1);
for i=1:nST
    ringFits{i}=single(0);
end;

iter=1;
%%
while iter<=nIters
    iteration=iter
    ri.useAltImgs=(iter==nIters) && useAltImgs;  % turn on use of alt images
    % at the last iteration
    %%
    allProjs=zeros(n0,n0,nRefs,2,nVols,nTwins,'single'); % 2 is sums/norms
    
    %    multiply by the volume mask
    refVols=refVols.*repmat(ri.volMask,1,1,1,nVols,nTwins);
    %     Set up twin volumes if necessary
    if fscFlag==1  % conventional fsc
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
    
    
    ri.useAltImgs=(iter==nIters) && useAltImgs;  % turn on use of alt images
    % at the last iteration
    %%
    allProjs=zeros(n0,n0,nRefs,2,nVols,nTwins,'single'); % 2 is sums/norms
    
    %    multiply by the volume mask
    refVols=refVols.*repmat(ri.volMask,1,1,1,nVols,nTwins);
    %     Set up twin volumes if necessary
    if fscFlag==1  % conventional fsc
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
                =reEMStep22(img1-ringFits{is},refs(:,:,:,:,iTw),si1{is},ri,mo1,altImg1);
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
                =reEMStep22(img1-ringFits{is},refs(:,:,:,:,iTw),si1{is},ri,mo1,altImg1);
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
    roi.pTransUnrot=zeros(nt,nt,nImgs,'single');
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
        roi.pTransUnrot(:,:,inds)=ro1{is}.pTransUnrot;
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
end;
    
