% re2DClassify3.m
% Distributed EM 2D classification
% Uses file-sharing to synchronize parallel jobs.
% The ri.timeout variable sets waiting times.
%
% We use three structures in the program.
% si "Stack info" is the structure, created by StackExtractor, describing
% the stack data including ctfs, active flags and original mi structures.
% ri "Run info" contains all global parameters, angle list, etc. for the
% reconstruction problem.
% moi "Model info" contains all model variables (incl. volumes)
% roi "Run output info" contains output from the E-Step including the
%       latent probabilities and class means.
% The collection of files from parallel computations occurs at two stages.
% - "group" roi structures (gRoi) are collected before reconstruction
% - "group" Fourier volumes (gfVs) are collected after the insertion step.
% The collected roi and reconstructed volumes are used to generate an moi
% file, which is then read (and the relevant group extracted) by each job.
doSimulateBatch=0;
allowRestart=1;
fixNGroups=1;
fixTwin=0;
frcExponent=.5;
nExponent=0;

% Figure out who we are
host=getenv('HOSTNAME');
isCluster=strncmp(host,'compute-',8);  % louise nodes start with 'compute-'

s=getenv('PBS_ENVIRONMENT');
sb=getenv('BATCH_ENVIRONMENT');
isInteractive=~(strcmp(s,'PBS_BATCH') || strcmp(sb,'BATCH') || doSimulateBatch);
showGraphics=isInteractive || doSimulateBatch;  % show graphics anyway if simulating.

numJobs=str2double(getenv('NUM_JOBS'));
if isnan(numJobs)
    numJobs=1;
end;
jobIndex=str2double(getenv('PBS_ARRAYID')); % zero-based index
if isnan(jobIndex)
    jobIndex=0;
end;

if doSimulateBatch
    numJobs=MyInput('numJobs',numJobs);
    jobIndex=MyInput('jobIndex',jobIndex);
end;


jobIndex
numJobs

if isCluster
    cd('/fastscratch/fjs2/151117/KvLipo80slot3');
else
%     We assume we're in the correct directory already.
%     cd('/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1')
    %     cd('/Users/fred/EMWork/Hideki/150114/KvLipo21slot1');
    %     cd('/Users/fred/EMWork/Hideki/150119/KvLipo14slot3');
end;


disp('Make run info');
[ri,si,moi]=re2DMakeRunInfo;  % --------construct ri, si, initial moi-----------
ri=re2DMakeRunInfoScaled(ri,ri.nCurrent);

%%
CheckAndMakeDir(ri.outPath,1);
CheckAndMakeDir(ri.tempPath,1);
CheckAndMakeDir([ri.outPath 'mrc/'],1);  % this will have class means
CheckAndMakeDir([ri.outPath 'jpeg/'],1);

if jobIndex==0
    riName=[ri.outPath 'ri.mat'];
    save(riName,'ri');
end;

% Initialize figures
if showGraphics
    figure(1);
    set(gcf,'position',[100,100,400,360]);
    figure(2);
    set(gcf,'position',[500,100,800,720]);
    figure(3);
    set(gcf,'position',[400,200,800,720]);
    figure(4);
end;

% in batch mode, numJobs is defined and is the total number of jobs used
% in this run, and jobIndex is in the range 0..numJobs-1
% iTwin is 1..ri.nTwins, i.e. 1..2
% To identify the current job, we have the variables
% nGroups, the number of jobs for each twin, and
% iGroup in the range 1..nGroups
if isInteractive
    if ri.nTwins>1 || fixTwin
        iTwin=MyInput('Twin value',1);
    else
        iTwin=1;
    end;
        ri.nGroups=MyInput('nGroups',ri.nGroups);
    nGroups=ri.nGroups;
    if ri.nGroups>1
        iGroup=MyInput('iGroup',1);
    else
        iGroup=1;
    end;
else  % batch mode.  numNodes must be even or else we have a duplicate group.
    iTwin=mod(jobIndex,ri.nTwins)+1;  % will be 1 or 2
    iGroup=floor(jobIndex/ri.nTwins)+1; % will be 1, 2, etc.
    nGroups=max(1,floor(numJobs/ri.nTwins));
    ri.nGroups=nGroups;
end;

% Create a log file.
logName=[ri.outPath char(iTwin+96) sprintf('%02d',iGroup) 'log.txt'];
logFile=fopen(logName,'a');  % Append to existing
logs=struct;
logs.handles=[1 logFile];
logs.idString='';
mdisp(logs,'--------------------------------------------------');
mprintf(logs.handles,'%s%s%d\n',datestr(now),' Startup: jobIndex ',jobIndex);
mprintf(logs.handles,'%s%s%s%s\n',' jobID ', jobID,' hostname ',host );
mprintf(logs.handles,'%s=%d  %s=%d  %s=%d\n',...
    'iTwin',iTwin,'iGroup',iGroup,'nGroups',ri.nGroups);
mprintf(logs.handles,'%s: %s\n','Log file',logName);
mdisp(logs,'Working dir', pwd);
mdisp(logs,'Output path',ri.outPath);

nTwin=ri.nTwin(iTwin);

% Put in the extensive model fields
n=size(moi.refs,1);
moi.imgAmps=moi.imgAmp*ones(nTwin,1,'single');
moi.ringFits=zeros(n,n,nTwin,'single');
moi.activeTrans=true(ri.nTrans^2,nTwin);
moi.activeRVs=true(ri.angleN(3),ri.angleN(2),ri.nVols);  % invariant

% Are we restarting?
[moiName,startIter]=reFindLatestMoi(ri.outPath,ri,iTwin);  % startIter=0 if nothing found.
if startIter>ri.startIter && allowRestart
    mdisp(logs,'Restarting at ',startIter);
    moi=LoadStruct([ri.outPath moiName]);
    ri=LoadStruct([ri.outPath 'ri.mat']); %%%%
    n=size(moi.refs,1);
    ri=re2DMakeRunInfoScaled(ri,n);
    ri.startIter=startIter;
    ri.nGroups=nGroups;
    %     Check for an roi file
    roiName=moiName;
    roiName(end-6)='r';  % moi.mat -> roi.mat
    if exist([ri.outPath roiName],'file')
        mdisp(logs,'Loading roi file, skipping EM step.');
        skipEM=1;
        roi=LoadStruct([ri.outPath roiName]);
    elseif ri.flags.restartRoiCollection % check to see if temp files exist
        if reGatherRois(ri,startIter,iTwin)
            mdisp(logs,'Gathering roi files, skipping EM');
            roi=reGatherRois(ri,iter,iTwin,logs,iGroup==1); % write only if iGroup=1
            skipEM=1;
        end;
    end;
else
    ri=reMakeRunInfoScaled(ri);
    if iGroup==1 % store our first moi
        save([ri.outPath reGetNameCode(ri,ri.startIter,iTwin) 'moi.mat'],'moi');
    end;
end;
frcTrend=[];
gImgs=0;  % start out with no images.
iter=ri.startIter;
%%  ---------------Iterations-----------------------------------------------
while iter<=ri.nIters
    %%    %     Have we already done an EM step?
    %     roiName=[reGetNameCode(ri,iter,iTwin) 'roi.mat'];
    %     [roi,ok]=reCheckAndLoadMat(ri.outPath, roiName,1);
    %     if ~ok  % roi doesn't already exist.
    
    %         Check and rescale the ri parameters
    iseq=find(ri.nSequence(:,2)>=iter,1);
    if numel(iseq)<1
        iseq=1;
    end;
    n=ri.nSequence(iseq,1);
        ri=re2DMakeRunInfoScaled(ri,n);

    if iGroup==1
        mdisp(logs,['Writing ' ri.outPath 'ri.mat']);
        save([ri.outPath 'ri.mat'],'ri');  % Save the latest ri structure
        if iter==ri.startIter  % if we're starting or restarting
            save([ri.outPath 'ri' sprintf('%02d',iter) '.mat'],'ri');
        end;
        
    mdisp(logs,'');
    mdisp(logs,'iter',iter);
    mdisp(logs,'Image size',n);
    
    %     Get the previous model
    if iter>ri.startIter && iGroup>1  % read moi from a file
        mdisp(logs,[datestr(now) ' Waiting for the moi file ']);
        [moi,ok]=reCheckAndLoadMat(ri.outPath, [reGetNameCode(ri,iter,iTwin) 'moi.mat'],...
            ri.timeout(1),1,numel(fieldnames(moi)));
        if ok
            mdisp(logs,[datestr(now) ' --found.']);
        else
            mdisp(logs,[datestr(now) ' --not found.']);
            fclose(logFile);
            return
        end;
    end;
    
    oldN=size(moi.refVols,1);
    if n~=oldN  % change of moi size needed
        moi=re2DRescaleMoi(moi,ri);
    end;
    if any(iter==ri.resetActiveRVIters)
        moi.activeRVs=true(size(moi.activeRVs));
        mdisp(logs,'Resetting moi.activeRVs');
    end;

    
    %     Normalize the ref amplitudes to unity variance
    refs=moi.refs;
    rSigma=sqrt(refs(:)'*refs(:)/numel(refs));
    moi.refs=refs/rSigma;
    
    %
    gFlags=reGetGroupFlags(ri,iTwin,iGroup);
    gMoi=reSplitStructureFields(moi,gFlags);
    
    %     Check for image-size consistency
    n=size(gImgs,1);  % the image size
    if iter==ri.startIter || ri.nCurrent ~= n
        % Reload the group's images
        loadAltImgs=ri.nCurrent==ri.nFinal; % are we at the final iterations' size?
        mdisp(logs,'Loading images.  Total images ', num2str(ri.nTwin(iTwin)));
        [gSi,gImgs,gAltImgs,groupSize]=reLoadStackGroup(ri,si,iTwin,iGroup,loadAltImgs);
        if gSi.pixA ~= ri.pixA
            warning(['Inconsistent pixel sizes ' num2str([gSi.pixA ri.pixA])]);
        end;
        n=size(gImgs,1);
        ngImgs=size(gImgs,3);
        mdisp(logs,'ngImgs loaded ',ngImgs);
    end;
    %
    % do the ring-fit correction
    imgs=gImgs-gMoi.ringFits;
    ri.usePrior= iter<3;  %%%%%%%%
    
    %%  ====================Do the E-step=========================
    mdisp(logs,[datestr(now) ' EM Step']);
    if iter<ri.nIters || ~ri.flags.useAltImgs  % only the last iter requires alt images
        [gRoi,refs]=reNodeEMStep(ri,gSi,gMoi,imgs,0,logs);
    else
        [gRoi,refs]=reNodeEMStep(ri,gSi,gMoi,imgs,gAltImgs-gMoi.ringFits);
    end;
    
    mdisp(logs,'inBoundImgs',sum(gRoi.inBoundImgs));
    
    %       Store the group roi structure
    if iGroup>1 || ri.flags.writeAllGroupFiles  % Encode iter, iGroup and iTwin in the output file name
        roiName=[ri.tempPath reGetNameCode(ri,iter,iTwin,iGroup) 'roi.mat'];
        save([roiName '_tmp'],'gRoi');
        eval(['!mv ' roiName '_tmp ' roiName]);
    end;
    %
    if showGraphics
        %     Display the group's latent probabilities
        figure(1);
        set(gcf,'position',[100,100,400,360]);
        figure(2);
        set(gcf,'position',[500,100,800,720]);
        reShowLatentVars(imgs,refs,ri,gMoi,gRoi,iter,[1 2]);
        drawnow;
        if ri.flags.saveFigures
            set(gcf,'paperpositionmode','auto');
            figName=[reGetNameCode(ri,iter,iTwin,iGroup) '.jpg'];
            print([ri.outPath 'jpeg/' figName],'-djpeg');
        end;
    end;
    
    %     Gather the roi structures
    if iGroup==1  % We are doing the master set
        roi=reGatherRois(ri,iter,iTwin,gRoi,logs);
        
        oldMoi=moi;
        
        if ri.nTwins>1  % get the other roi
            xRoiName=[reGetNameCode(ri,iter,3-iTwin) 'roi.mat'];
            mdisp(logs,[datestr(now) ' Loading the other roi file ' xRoiName]);
            [xRoi,ok]=reCheckAndLoadMat(ri.outPath, xRoiName,ri.timeout(2),...
                1,numel(fieldnames(roi)));
            
            %        Compute the frcs
            disp('Computing FRCs');
            nCls=size(roi.classMeans,3);
            frc=FRC(roi.classMeans,xRoi.classMeans);
            h=reFitErf(frc);
            h=max(0,(h-h(end))/(h(1)-h(end))).^frcExponent;  % normalize the filter
            hRot=ToRect(h);
            k0=1e-3;
            k1=n^nExponent;
            kx=k0+(1-hRot)./(hRot+k0);
            kx=repmat(kx,1,1,nCls)*k1;
            roi.classMeans=roi.classMeans+xRoi.classMeans;
            roi.classNorms=roi.classNorms+xRoi.classNorms;
            frcTrend(1:n/2,iter)=frc;
            figure(4);
            fs=(0:n/2-1)/(ri.pixA*n);
            plot(fs,frcTrend,'-',fs,frc,'k-',fs,h,'k--');
        else
            kx=n^nExponent;
            disp(['Using constant k' num2str(kx)]);
        end;
        
        %%     =====================Do the 2D reconstruction======================
        
        fCls=fftshift2(fft2(roi.classMeans));
        nCls=size(fCls,3);
%         %     Wiener normalization
%         f=RadiusNorm(n)/ri.pixA;
%         f0=.2/ri.pixA;
%         k=.001*(1+exp((f/f0).^2));
%         kx=repmat(k,1,1,nCls);
%         %         kx=.001;
        moi.refs=real(ifft2(ifftshift2((fCls)./(kx*moi.sigmaN^2+roi.classNorms))));
        sds=sqrt(mean(mean(moi.refs.^2)));
        k0=.001;
        moi.refs=moi.refs./repmat(sds+k0,n,n,1);
        %
        %
        % Assign other model values
        nImgs=numel(roi.imgAmps);
        moi.sigmaN=sqrt(mean(roi.varNs));
        moi.sigmaC=oldMoi.sigmaC;
        moi.sigmaG=oldMoi.sigmaG;
        minPRef=1e-4;
        moi.pRefs=max(roi.pRefs,minPRef);
        moi.imgAmps=median(abs(roi.imgAmps))*ones(nImgs,1);
        %         reSetImageAmplitudes2(si,ri,iTwin,roi.imgAmps);
        moi.ringFits=oldMoi.ringFits+roi.ringFits;
        moi.activeTrans=reAssignActiveTrans(roi.pTrans,ri.thresholds(1));
        moi.pVols=mean(roi.pVols,2);
        moi.intensiveFields={'sigmaN' 'sigmaC' 'sigmaG' 'refVols', 'pVols','pRefs','refs'};
        
        fracActiveTrans=sum(moi.activeTrans(:))/numel(moi.activeTrans);
        mdisp(logs,'fracActiveTrans',fracActiveTrans);
        %     fracActiveAlphas=sum(roi.pAlphas(:)>transThreshold)/numel(roi.pAlphas)
        %     fracActiveRefs=sum(roi.pRefs(:)>transThreshold)/numel(roi.pRefs)
        moiName=[ri.outPath reGetNameCode(ri,iter+1,iTwin) 'moi.mat'];
        save([moiName '_tmp'],'moi');
        eval(['!mv ' moiName '_tmp ' moiName]);
        
        %     Display the volumes
        if showGraphics
            figure(3);
            %             set(gcf,'position',[100,300,400,400]);
            ImagicDisplay2(moi.refs);
            title(figName);
            set(gcf,'paperpositionmode','auto');
            print([ri.outPath 'jpeg/' figName],'-djpeg');
        end;
        
        
        
    end;      % if iGroup==1
    iter=iter+1;
end;          % while iter
fclose(logFile);
