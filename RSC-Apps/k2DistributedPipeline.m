% k2DistributedPipeline.m
% Run the processing pipeline for K2 or F2 movies.
% This is assumed to run after running k2CreateInfoFiles or f2CreateInfoFiles
% on one worker so
% we only read the mi files, not create them.

f2Mode        =0;  % 0 means k2 data
serialMode    =0;   % go through all the steps before moving to next micrograph
% -- currently doesn't work with serialMode=0 !
%checkLogs     =0;  % not yet in use.

findUnfinished=0; % Don't process data, just make a new allNamesUnf.mat file.
% This option makes sense only if nJobs=1, serialMode is used.

doFindJump    =0;
doTrack       =0;  % do movie alignment
doMerge       =0;
forceMerging  =0;
doDownsampleMergedImages=0;
doCompressMovies      =0;  % compress movies
doCompressMicrographs =0;  % compress micrographs
doFindVesicles        =0;
%*** special multiple VesicleFinder runs ****
doMultiFindVesicles   = 0;
% findVesicleAmps=[5e-4 6e-4 7e-4 8e-4];
% findVesicleDirs={'KvLipo121_2w10_v3.5/';'KvLipo121_2w10_v3.6/';
%                  'KvLipo121_2w10_v3.7/';'KvLipo121_2w10_v3.8/'};
%***
doPrelimInverseFilter =0;
doRefineVesicles      =1;
forceRefineVesicles   =0;
refineVesicleAmpsOnly=0;
minRefineVesiclesSequence=0;  % 0 if don't consider.
% minRefineVesiclesSequence=inf ;    % inf forces refinement
doInverseFilter       =0;
forceInverseFilter=0;
maxAge=1;  % Re-run if the old file is older than this number of days.
% days before the present we go ahead and re-run the
% function.  So, to re-run processing if the latest log entry is < 1 day old,
% set minAge=1.

doPickingPreprocessor =1;


workingDir='/gpfs/ysm/scratch60/sigworth/fjs2/200707/';
workingDir='/gpfs/ysm/scratch60/sigworth/hs468/DataFromRIKEN/200816/025015_1_1/';
workingDir='/gpfs/ysm/scratch60/sigworth/fjs2/20211122/';
infoDir='Info_xchg/';

compressedDir=[workingDir 'Compressed/'];
localWorkingDir=workingDir;
logDir='Log/';

createMiFiles=0;

pars=struct;
pars.overwrite=0;
pars.useParfor=0;

pars.forceTracking=0;  % Drift Tracker
pars.writeZTiff=0;     % 2= z.tif only; 1= write both .mrc and z.tif
pars.writeGraphics=1;
pars.testSegments=[2 20; 38 inf];
pars.testSegments=[2 20; 25 inf]; % 161101 data
pars.writeStack=1;  % Dirft tracker

% for merging
%pars.defaultPixA=1.05;
%pars.defaultPixA=1.781;
%pars.defaultPixA=0;
pars.searchDefoci=[1 8 ; 8 15]; % for MergeImages [1stmin 2ndMin ; 1stMax 2ndMax]
pars.doAlignment=0;  % MergeImages
pars.doFitting=0;
pars.doWriteInfos=1;
%pars.weights=[1 0];  %%% single exposure
%pars.weights=[1 1];
pars.weights=1;
pars.mcDS=1;
pars.mergeMode=3;  %%% no phase flip!
%pars.mergeMode=1;   %%% normal
pars.mapMode='Kv';

pars.UsePWFilter=doPrelimInverseFilter;

% ---for rsRefineVesicleFits----
pars.doPreSubtraction=1;  % rsRefineVesicleFits: pre-subtract old vesicle fit.
% pars.rTerms=[100 150 200 300  inf];
pars.rTerms=[90 100 120 150 200 250 300 inf];
pars.dsSmall=4; % downsampling of 'small' merged image
pars.overwrite=0;

pars.loadFilenames=0; % pick up allNames.mat in base directory
pars.cpe=0;  % 0 means no change.

% pars.modelMiName='~/scratch60/170417/KvLipo134_4/sq02w11/Info/sq02_1_0001_Apr18_15.11.43mi.txt';
%pars.modelMiName='~/scratch60/170609/KvLipo135_1/Info/sq02_1_0001_Jun09_17.25.06mi.txt';
% pars.modelMiName='/ysm-gpfs/scratch60/fjs2/170814/KvLipo125_3a/Info/sq05_1_0001_Aug14_14.58.28mi.txt';
pars.modelMiName='~/data/MembraneRef/160909_sq02_1_01mi.txt';


doSimulateBatch=0;  % for simulating batch on local machine

% Figure out who we are
host=getenv('HOSTNAME');
isCluster=strncmpi(host,'c',1);  % farnam nodes are like 'c22n09'
disp(['isCluster = ' num2str(isCluster)]);
sb=getenv('BATCH_ENVIRONMENT');
isInteractive=~( strcmp(sb,'BATCH') || doSimulateBatch );
disp(['isInteractive = ' num2str(isInteractive)]);
pars.showGraphics=isInteractive || doSimulateBatch;  % show graphics anyway if simulating.
pars.showGraphics=0; %%%%%%%%

numJobs=str2double(getenv('NUM_JOBS'));
jobIndex=str2double(getenv('JOB_ID')); % one-based index

% Set the current datstring
refDate=now;

% %%%%%%%%%%%%%%%%%%%%%
% if jobIndex==12  % special case!!!
% doTrack       =1;  % do movie alignment
% doMerge       =1;
% end;
% %%%%%%%%%%%%%%%%%%%%%%%%

if doSimulateBatch || isnan(numJobs) || isnan(jobIndex)
    if isnan(numJobs)
        numJobs=1;
    end;
    if isnan(jobIndex)
        jobIndex=1;
    end;
    
    disp('Batch simulation:');
    numJobs=MyInput('numJobs',numJobs);
    jobIndex=MyInput('jobIndex',jobIndex);
end;

if isCluster
    cd(workingDir);
else
    cd(localWorkingDir)
end;

% Set up the log file
CheckAndMakeDir('Log',1);

% Create a log file.
logName=[logDir sprintf('%02d',jobIndex) 'log.txt'];
logFile=fopen(logName,'a');
logs=struct;
logs.handles=[1 logFile];
logs.idString='';
pars.logs=logs;
mdisp(pars.logs,' ');
mdisp(pars.logs,'======================================================================-');
mprintf(pars.logs.handles,'%s%s%d\n',datestr(now),' Startup: group ',jobIndex);
mdisp(pars.logs,pwd);

% Get the mi file names
% %if f2Mode
%     allNames=f2FindInfoFiles;
%     nNames=numel(allNames);
% % else
% %     pars0=pars;
% %     pars0.overwrite=0;  % don't replace extant files.
% %     allNames=k2CreateInfoFiles(pars0);
% % end;
if pars.loadFilenames
    disp('loading filenames');
    load allNames.mat
else
    disp('Finding the mi files');
    allNames=f2FindInfoFiles(infoDir);
end;
nNames=numel(allNames);
disp([num2str(nNames) ' files total']);

if nNames<1
    msg=['No mi files found in ' pwd '/Info/'];
    mdisp(pars.logs,msg);
    error(msg);
end;

ngx=0; % extra names to worker 1

% Figure out what we will do
blockSize=(nNames)/numJobs;
% group 1 has an extra ngx entries.
ourBlockStart=round((jobIndex-1)*blockSize)+1+ngx*(jobIndex>1);
ourBlockEnd=min(nNames,round(jobIndex*blockSize)+ngx);
mprintf(pars.logs.handles,'files %d to %d\n',ourBlockStart,ourBlockEnd);

allNamesUnf=cell(0,1); % unfinished files
nUnfinished=0;

jobNames=allNames(ourBlockStart:ourBlockEnd);
numJobNames=numel(jobNames);
if numJobNames<1  % nothing to do
    error(['No mi files found: ' pwd]);
end;

if serialMode
    nNames=1;
else
    nNames=numJobNames;
    findUnfinished=0; % can do this only in serial Mode
end;
iName=(1:nNames);
while iName(end)<=numJobNames
    disp(['Working on images ' num2str(iName(1)) ' to ' num2str(iName(end)) ...
        ' of ' num2str(numJobNames)]);
    unfinished=0;
    
    ourNames=jobNames(iName);
    if serialMode
        disp(ourNames{1});
        mi=ReadMiFile(ourNames{1});
        [logSequence,dates]=miDecodeLog(mi);
    else
        logSequence=true(1,10);
        dates=zeros(1,10);
    end;
    
    % find jump (sequence 1)
    if doFindJump && ~f2Mode
        if findUnfinished
            unfinished=1 | unfinished;
        else
            k2FindDefocusJump(ourNames,pars);
        end;
    end;
    
    % drift tracker (sequence 2)
    if doTrack
        if findUnfinished
            unfinished=1 | unfinished;
        else
            if f2Mode
                f2DriftTracker(ourNames,pars);
            else
                k2DriftTracker(ourNames,pars);
            end;
        end;
    end;
    % merge images (sequence 3)
    if doMerge
        if ~logSequence(3) || forceMerging || now-dates(3)>maxAge
            if findUnfinished
                unfinished=1 | unfinished;
            else
                MergeImages(ourNames,pars);
            end;
        else
            disp('Merging skipped.');
        end;
    end;
    
    if ~ findUnfinished
        doCompressMovies=doCompressMovies && f2Mode;  % can't do this with k2 movies.
        if doCompressMovies || doCompressMicrographs
            f2CompressMovies(ourNames,compressedDir,doCompressMovies,doCompressMicrographs);
        end;
        if doDownsampleMergedImages
            DownsampleMergedImages(ourNames);
        end;
    end;
    % inverse filter (sequence 6)
    if doPrelimInverseFilter && logSequence(6)<=logSequence(3)
        fpars=struct;
        fpars.useUnsubImage=1;
        if findUnfinished
            unfinished=1 | unfinished;
        else
            meInverseFilterAuto(ourNames,fpars);
        end;
    end;
    % find vesicles (sequence 4) *****************
    if doFindVesicles && (logSequence(4)<=logSequence(3) || now-dates(4)>maxAge)
        if findUnfinished
            unfinished=1 | unfinished;
        else
            %            VesicleFinder(ourNames);
            Vesicle_finding_GUI(ourNames);
        end;
    elseif doMultiFindVesicles && ~findUnfinished
        for i=1:numel(findVesicleAmps)
            vfpars.sav.vesicleAmps=[findVesicleAmps(i) 2e-3 0];
            cd('..');
            cd(findVesicleDirs{i});
            disp(pwd);
            VesicleFinder(ourNames,vfpars);
        end;
    end;
    % refine vesicles (sequence 5)
    if forceRefineVesicles ...
            || (doRefineVesicles && (logSequence(5)<logSequence(4) ...
            || logSequence (5)< minRefineVesiclesSequence || now-dates(5)>maxAge))
        rpars=pars;
        if findUnfinished
            unfinished=1 | unfinished;
        else
            if refineVesicleAmpsOnly
                rpars.fitModes={'LinOnly'};
                rpars.fractionStartingTerms=1; % total terms to use in each round
                rpars.fractionAmpTerms=1;
                % Extra peaks in the scattering profile
                rpars.peakPositionA=[-37 0 37];  % empirical default.  Works a bit better than [37 37]
                rpars.targetPixA=10;  % downsampled image resolution for radius fitting
                
                rpars.xPeakSigmaA={5 5}; % width of extra Gaussian peaks, in angstrom
                %     The following must have at least as many elements as dpars.fitModes!
            end;
            
            rsRefineVesicleFits(ourNames,rpars);
            if serialMode  % update the log sequence
                mi=ReadMiFile(ourNames{1});
                [logSequence,dates]=miDecodeLog(mi);
            end;
        end;
    elseif doRefineVesicles && ~findUnfinished
        disp('  Refine Vesicles skipped.');
    end;
    % inverse filter (sequence 6)
    if doInverseFilter && (logSequence(6) <= logSequence(5) ...
            || now-dates(6)>maxAge)
        if findUnfinished
            unfinished=1 | unfinished;
        else
            meInverseFilterAuto(ourNames);
        end;
    elseif doInverseFilter
        disp('  Inverse Filter skipped.');
    end;
    % picking preprocessor (sequence 8)
    if doPickingPreprocessor && (logSequence(8) <= logSequence(6) || ...
            now-dates(8)>maxAge)
        if findUnfinished
            unfinished=1 | unfinished;
        else
            % no picking after latest vesicle refinement? Then run it.
            rsPickingPreprocessor4(ourNames,pars);
        end;
    end;
    
    if findUnfinished && unfinished
        nUnfinished=nUnfinished+1;
        disp(['-----unfinished ' num2str(nUnfinished)]);
        allNamesUnf(nUnfinished,1)=ourNames;
    end;
    iName=iName+nNames;
end;

if findUnfinished
    disp(['Saving ' num2str(nUnfinished) ' unfinished file names']);
    allNames=allNamesUnf;
    save allNamesUnf.mat allNames ;
end;

