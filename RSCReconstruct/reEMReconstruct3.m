% reEMReconstruct3.m
% Distributed EM reconstruction
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

% Figure out who we are
host=getenv('HOSTNAME');
isCluster=strncmp(host,'compute-',8);  % louise nodes start with 'compute-'

s=getenv('PBS_ENVIRONMENT');
sb=getenv('BATCH_ENVIRONMENT');
isInteractive=~(strcmp(s,'PBS_BATCH') || strcmp(sb,'BATCH') || doSimulateBatch);
% isInteractive=~(strcmp(s,'PBS_BATCH') || doSimulateBatch);
showGraphics=isInteractive || doSimulateBatch;  % show graphics anyway if simulating.

numJobs=str2double(getenv('NUM_JOBS'));
if isnan(numJobs)
    numJobs=1;
end;
jobIndex=str2double(getenv('PBS_ARRAYID')); % zero-based index
if isnan(jobIndex)
    jobIndex=0;
end;
jobID=getenv('PBS_JOBID');

if doSimulateBatch
    numJobs=MyInput('numJobs',numJobs);
    jobIndex=MyInput('jobIndex',jobIndex);
end;


jobIndex
numJobs

if isCluster
    basePath='/fastscratch/fjs2/';
else
    basePath='/Users/fred/EMWork/Hideki/';
end;

%cd([basePath '140625n');
cd([basePath '151117/KvLipo80slot3/']);
%cd([basePath '160622Sim/']);
%
disp('Make run info');
[ri,si,moi]=reMakeRunInfo;  % --------construct ri, si, initial moi-----------


%%%%%special restart-------------------------------------------------------
%load('Reconstructions/Recon112n1/ri.mat'); %%%%%%%%%%%%%%%%%%
% load('Reconstructions/Recon112ct/i11a_moi.mat');

% moi.imgAmp=median(moi.imgAmps);
% figure(3);
% ShowSections(moi.refVols,[],45);
% ri.nCurrent=48;
% ri.nIters=12;
% ri.activeFlags(1000:end)=false;
% ri.twinFlags=[ri.activeFlags ri.activeFlags];
% ri.nTwin=sum(ri.twinFlags);
% %%%%%%%%%

ri=reMakeRunInfoScaled(ri,ri.nCurrent);


CheckAndMakeDir(ri.outPath,1);
CheckAndMakeDir(ri.tempPath,1);
CheckAndMakeDir([ri.outPath 'mrc/'],1);
CheckAndMakeDir([ri.outPath 'jpeg/'],1);

if jobIndex==0 && ~ri.flags.restartRoiCollection
    str=['!rm ' ri.tempPath '*'];
    disp(str);
    eval(str);
end;
%     riName=[ri.outPath 'ri.mat'];
%     save(riName,'ri');
% end;
%
% Initialize figures
if showGraphics
    figure(1);
    %     set(gcf,'position',[100,100,400,360]);
    figure(2);
    %     set(gcf,'position',[500,100,800,720]);
    figure(3);
    %     set(gcf,'position',[400,200,800,720]);
end;

% in batch mode, numJobs is defined and is the total number of jobs used
% in this run, and jobIndex is in the range 0..numJobs-1
% iTwin is 1..ri.nTwins, i.e. 1..2
% To identify the current job, we have the variables
% nGroups, the number of jobs for each twin, and
% iGroup in the range 1..nGroups
if isInteractive
    if ri.nTwins>1
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
n=size(moi.refVols,1);
moi.imgAmps=moi.imgAmp*ones(nTwin,1,'single');
moi.activeTrans=true(ri.nTrans^2,nTwin);
moi.activeAlphas=true([size(ri.alphas) nTwin]);
moi.activeRVs=true(ri.angleN(3),ri.angleN(2),ri.nVols,nTwin);

% Are we restarting?
skipEM=0;
[moiName,startIter]=reFindLatestMoi(ri.outPath,ri,iTwin);  % startIter=0 if nothing found.
if startIter>ri.startIter && ri.flags.allowRestart
    mdisp(logs,'Restarting at ',startIter);
    moi=LoadStruct([ri.outPath moiName]);
    %     ri=LoadStruct([ri.outPath 'ri.mat']); %%%%
    n=size(moi.refVols,1);
    ri=reMakeRunInfoScaled(ri,n);
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

gImgs=0;  % start out with no images.

%%  ---------------Iterations-----------------------------------------------
for iter=ri.startIter:ri.nIters
    %     Have we already done an EM step?
    %     roiName=[reGetNameCode(ri,iter,iTwin) 'roi.mat'];
    %     [roi,ok]=reCheckAndLoadMat(ri.outPath, roiName,1);
    %     if ~ok  % roi doesn't already exist.
    %         Find the new image size
    iseq=find(ri.nSequence(:,2)>=iter,1);
    if numel(iseq)<1
        iseq=1;
    end;
    n=ri.nSequence(iseq,1);  % image size for this iteration
    ri=reMakeRunInfoScaled(ri,n);
    %
    %         oldN=size(moi.refVols,1);
    %         if n~=oldN  % change of image size
    %             ri=reMakeRunInfoScaled(ri,n);
    %             moi=reRescaleMoi(moi,ri);
    %         end;
    if iGroup==1
        mdisp(logs,['Writing ' ri.outPath 'ri.mat']);
        save([ri.outPath 'ri.mat'],'ri');  % Save the latest ri structure
        if iter==ri.startIter  % if we're starting or restarting
            save([ri.outPath 'ri' sprintf('%02d',iter) '.mat'],'ri');
        end;
    end;
    
    mdisp(logs,'');
    mdisp(logs,'iter',iter);
    mdisp(logs,'Image size',n);
    %     Get the previous model if we're not starting up
    if iter>ri.startIter && iGroup>1  % read moi from a file
        mdisp(logs,[datestr(now) ' Waiting for the moi file ']);
        [moi,ok]=reCheckAndLoadMat(ri.outPath, [reGetNameCode(ri,iter,iTwin) 'moi.mat'],...
            ri.timeout(1),1,numel(fieldnames(moi)));
        if ok
            mdisp(logs,[datestr(now) ' --found.']);
        else
            mdisp(logs,[datestr(now) ' --not found. Halting.']);
            fclose(logFile);
            return
        end;
    end;
    oldN=size(moi.refVols,1);
    if n~=oldN  % change of moi size
        moi=reRescaleMoi(moi,ri);
    end;
    if any(iter==ri.resetActiveRVIters)
        moi.activeRVs=true(size(moi.activeRVs));
        mdisp(logs,'Resetting moi.activeRVs');
    end;
    
    if ~skipEM
        
        gFlags=reGetGroupFlags(ri,iTwin,iGroup);
        gMoi=reSplitStructureFields(moi,gFlags);
        
        %     Check for image-size consistency
        n=size(gImgs,1);  % the image size
        if iter==ri.startIter || ri.nCurrent ~= n % || iter>=ri.altImgIter
            % Load the group's images
            loadAltImgs=iter>=ri.altImgIter; % are we at the final iterations?
            mdisp(logs,'Loading images.  Total images ', num2str(ri.nTwin(iTwin)));
            [gSi,gImgs,gAltImgs,groupSize]=reLoadStackGroup(ri,si,iTwin,iGroup,loadAltImgs);
            if gSi.pixA ~= ri.pixA
                warning(['Inconsistent pixel sizes ' num2str([gSi.pixA ri.pixA])]);
            end;
            n=size(gImgs,1);
            ngImgs=size(gImgs,3);
            mdisp(logs,'ngImgs loaded ',ngImgs);
        end;
        
        %  ====================Do the E-step=========================
        mdisp(logs,[datestr(now) ' EM Step ' num2str(size(gImgs,3)) ' images']);
        if iter<ri.altImgIter || ~ri.flags.useAltImgs  % only the last iter requires alt images
            [gRoi,refs]=reNodeEMStep(ri,gSi,gMoi,gImgs,0,logs);
        else
            mdisp(logs,[datestr(now) ' using alt images ' num2str(size(gAltImgs,3))]);
            %            [gRoi,refs]=reNodeEMStep(ri,gSi,gMoi,imgs,gAltImgs-gMoi.ringFits);
            [gRoi,refs]=reNodeEMStep(ri,gSi,gMoi,gImgs,gAltImgs,logs);
        end;
        
        mdisp(logs,'inBoundImgs',sum(gRoi.inBoundImgs));
        
        %       Store the group roi structure
        roiName=[ri.tempPath reGetNameCode(ri,iter,iTwin,iGroup) 'roi.mat'];
        save([roiName '_tmp'],'gRoi');
        eval(['!mv ' roiName '_tmp ' roiName]);
        
        %%
        if showGraphics
            %     Display the group's latent probabilities
            figure(1);
            set(gcf,'position',[100,100,400,360]);
            figure(2);
            set(gcf,'position',[500,100,800,720]);
            reShowLatentVars(gImgs,refs,ri,gMoi,gRoi,iter,[1 2]);
            drawnow;
            if ri.flags.saveFigures
                set(gcf,'paperpositionmode','auto');
                figName=[reGetNameCode(ri,iter,iTwin,iGroup) '.jpg'];
                print([ri.outPath 'jpeg/' figName],'-djpeg');
            end;
        end;
        
        %     Gather the roi structures
        if iGroup==1  % We are doing the master set
            roi=reGatherRois(ri,iter,iTwin,logs);
        else
            roi=0;
        end;
    else
        skipEM=0;
    end;  % if roi already exists
    
    
    %%     =====================Do the reconstruction======================
    isStartup=iter<3;  % first two iterations we just re-est sigmaN and imgAmp
    
    if ri.flags.mode2D
        if iGroup==1  % do the reconstruction
            moi=reEstimateModelValues(ri,si,roi,iTwin,moi,isStartup,logs);  % updates everything except volumes.
            fCls=fftshift2(fft2(roi.classMeans));

            %         %     Wiener normalization
%         nCls=size(fCls,3);
%         f=RadiusNorm(n)/ri.pixA;
%         f0=.2/ri.pixA;
%         k=.001*(1+exp((f/f0).^2));
%         kx=repmat(k,1,1,nCls);
         kx=1;

            moi.refs=real(ifft2(ifftshift2((fCls)./(kx*moi.sigmaN^2+roi.classNorms))));
            sds=sqrt(mean(mean(moi.refs.^2)));
            k0=1;
            moi.refs=moi.refs./repmat(sds+k0,n,n,1);
            if showGraphics
                figure(10);
                imagsar(moi.refs);
            end;
        end;

    else  % 3D Reconstruction        
        if iGroup==1
            oldMoi=moi;
            moi=reEstimateModelValues(ri,si,roi,iTwin,moi,isStartup,logs);  % updates everything except volumes.
            mdisp(logs,['sigmaN=' num2str(moi.sigmaN) '  imgAmps=' num2str(mean(moi.imgAmps))]);
            vols=moi.refVols;  % default volumes
        end;
        if ~isStartup %  Do the reconstruction only after noise and amp estimates have settled.
            %  On every worker, do the distributed Fourier insertion
            gFvs=reGroupInsert(ri,iter,iTwin,iGroup,roi,logs);
            %
            if iGroup==1 % Do the overall reconstruction
                %   -----group 1:   Gather the Fourier insertions and normalize----
                [moi,vols]=reGroupReconstruct(ri,si,iter,iTwin,roi,moi,gFvs,logs);
                mdisp(logs,'Masking the references with ri.volMask');
                for i=1:ri.nVols
                    moi.refVols(:,:,:,i)=moi.refVols(:,:,:,i).*ri.volMask;
                end;
            end;
        end;
        
        if iGroup==1
            % Save the moi structure
            moiName=[ri.outPath reGetNameCode(ri,iter+1,iTwin) 'moi.mat'];
            save([moiName '_tmp'],'moi');
            system(['mv ' moiName '_tmp ' moiName]);
        end;
        
        if iGroup==1 && ~isStartup
            if showGraphics
                for iVol=1:ri.nVols
                    figure(iVol+2);
                    set(gcf,'position',[100,300*200*iVol,400,400]);
                    ShowSections(vols(:,:,:,iVol),[],45);
                    if ri.flags.saveFigures
                        figName=[reGetNameCode(ri,iter,iTwin,-iVol) '.jpg'];
                        subplot(3,3,1);
                        title(figName);
                        set(gcf,'paperpositionmode','auto');
                        print([ri.outPath 'jpeg/' figName],'-djpeg');
                    end;
                end;
            end;
            
            % Save our volumes
            for iVol=1:ri.nVols
                volName=[reGetNameCode(ri,iter,iTwin,-iVol) '.mrc'];
                WriteMRC(vols(:,:,:,iVol),ri.pixA,[ri.outPath 'mrc/' volName]);
            end;
            %%
            if ri.nTwins>1 %      Look up the other twin's volumes and compute the fsc
                moi.fscs=zeros(n/2,ri.nVols);
                otherVols=zeros(n,n,n,ri.nVols,'single');
                mdisp(logs,[datestr(now) ' Loading the other twin''s volumes.']);
                
                if iter>=ri.mergeIter
                    twinVolTimeout=ri.timeout(3);
                else
                    twinVolTimeout=ri.timeout(3)/10;
                end;
                
                for iVol=1:ri.nVols
                    altName=[reGetNameCode(ri,iter,3-iTwin,-iVol) '.mrc'];
                    [v,s,ok]=reCheckAndLoadMRC([ri.outPath 'mrc/'], altName,twinVolTimeout,5);
                    if ~ok
                        mdisp(logs,[datestr(now) ' load failed']);
                        break;
                    else
                        %                     if size(v,1)~=n % different voxel size
                        %                         v=DownsampleGeneral(v,n);
                        %                     end;
                        %                     otherVols(:,:,:,iVol)=v;
                        if all(size(otherVols)==size(vols))
                            
                            %                 Write the merged volume
                            volName=[reGetNameCode(ri,iter,3,-iVol) '.mrc'];  % volume 'c'
                            WriteMRC(otherVols(:,:,:,iVol)+vols(:,:,:,iVol),ri.pixA,[ri.outPath 'mrc/' volName]);
                            %                 Save the figure again, this time with the fsc
                            moi.fscs(:,iVol)=FSCorr2(otherVols(:,:,:,iVol),vols(:,:,:,iVol));
                            if showGraphics
                                figure(iVol+2)
                                subplot(3,3,9);
                                plot(moi.fscs(:,iVol));
                                if ri.flags.saveFigures
                                    figName=[reGetNameCode(ri,iter,iTwin,-iVol) '.jpg'];
                                    subplot(3,3,1);
                                    title(figName);
                                    set(gcf,'paperpositionmode','auto');
                                    print([ri.outPath 'jpeg/' figName],'-djpeg');
                                end;
                            end;
                        end;
                    end;
                end;  % for iVol
            else % nTwins=1, just save the moi file
                moiName=[ri.outPath reGetNameCode(ri,iter+1,iTwin) 'moi.mat'];
                save([moiName '_tmp'],'moi');
                system(['mv ' moiName '_tmp ' moiName]);
                
            end;
            if iter>=ri.mergeIter && ri.nTwins==2 % we'll merge the moi refVols
                if ~ok  % failure in reading fthe other volumes
                    %                 mdisp(logs,'Couldn''t create the twin-merged volumes.  .');
                    %                 fclose(logFile);
                    warning('Couldn''t create the twin-merged volumes, continuing anyway.');
                else
                    otherRefVols=reNormalizeModels(ri,otherVols);
                    if iTwin==1
                        moi.refVols=(reAlignVolumes(moi.refVols,otherRefVols)+moi.refVols)/2;
                    else
                        moi.refVols=(reAlignVolumes(otherRefVols,moi.refVols)+otherRefVols)/2;
                    end;
                    %                 moi.refVols=(moi.refVols+otherRefVols)/(2);
                    %             Save the moi again with the merged reference volumes
                    moiName=[ri.outPath reGetNameCode(ri,iter+1,iTwin) 'moi.mat'];
                    save([moiName '_tmp'],'moi');
                    system(['mv ' moiName '_tmp ' moiName]);
                end;
            end;
            if iter==ri.skewIter
                mdisp(logs,'Skew modification of moi.refVols');
                for i=1:ri.nVols
                    sls=round(.35*n):round(.42*n);
                    moi.refVols(:,:,:,i)=SkewVolume(moi.refVols(:,:,:,i),sls,30);
                    moiName=[ri.outPath reGetNameCode(ri,iter+1,iTwin) 'moi.mat'];
                    save([moiName '_tmp'],'moi');
                    system(['mv ' moiName '_tmp ' moiName]);
                end;
            end;
        end;  % if iGroup==1 && ~isStartup
    end;  % if mode2D
end;          % for iter

if iGroup==1  % Let's wait till all the twins are done
    theTwins=1:ri.nTwins;
else
    theTwins=iTwin;
end;
for iTw=theTwins
    mdisp(logs,[datestr(now) ' Waiting for the final moi file, twin= ' num2str(iTw)]);
    [moi,ok]=reCheckAndLoadMat(ri.outPath, [reGetNameCode(ri,iter,iTw) 'moi.mat'],...
        ri.timeout(1),1,numel(fieldnames(moi)));
    if ok
        str=' ...successful.';
    else
        str=' ...timed out.';
    end;
    mdisp(logs,[datestr(now) str]);
end;
fclose(logFile);
