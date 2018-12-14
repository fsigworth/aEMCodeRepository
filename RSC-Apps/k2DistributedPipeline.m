% k2DistributedPipeline.m
% Run the processing pipeline for K2 or F2 movies.
% This is assumed to run after running k2CreateInfoFiles or f2CreateInfoFiles
% on one worker so
% we only read the mi files, not create them.

f2Mode        =0;  % 0 means k2 data
serialMode    =1;   % go through all the steps before moving to next micrograph
%checkLogs     =0;  % not yet in use.

maxAge=1;         % Runs the operation is the previous one is older than this number of days.

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
doRefineVesicles      =0;
refineVesicleAmpsOnly=0;
%%% minRefineVesiclesSequence=0;  % 0 if don't consider.
minRefineVesiclesSequence=0;    % inf forces refinement
doInverseFilter       =1;
forceInverseFilter=1;
minAge=.01;  % if the corresponding log entry has a date stamp < minAge
% days before the present we go ahead and re-run the
% function.  So, to re-run processing if the latest log entry is < 1 day old,
% set minAge=1.

doPickingPreprocessor =0;


workingDir='/gpfs/ysm/project/fjs2/180226/Kv_1/';
%workingDir='/ysm-gpfs/pi/cryoem/krios/20171120/KvLipo123_1/'
%workingDir='/ysm-gpfs/scratch60/fjs2/160909/KvLipo121_2w11v3m1/'
% workingDir='/ysm-gpfs/scratch60/fjs2/160909/KvLipo121_2w10v3t/'
% workingDir='/ysm-gpfs/scratch60/fjs2/170808p/SimpleVes_raFit/';
%workingDir='~/project/180226/Kv_1selw10/';
%workingDir='~/project/20180620/';
%workingDir='~/project/20181025/20Frames/';
%workingDir='/gpfs/ysm/scratch60/fjs2/180226/Kv_1SelW10/';
%workingDir='/ysm-gpfs/scratch60/fjs2/170926Nelli/'
%workingDir='/ysm-gpfs/scratch60/fjs2/171031Nelli/'
%workingDir='/ysm-gpfs/scratch60/fjs2/170808/SiW10/'
%workingDir='/ysm-gpfs/scratch60/fjs2/170417/KvLipo134_4/sq04w11/';
%workingDir='/fastscratch/fjs2/160909/KvLipo121_2retw10ds1/'
%workingDir='/net/scratch2/fjs2/transfer/161124/KvW366FLipo_1/';
%workingDir='/net/scratch2/fjs2/transfer/161124/KvLipo/';
%workingDir='/Volumes/D215/160909/KvLipo121_2';
%workingDir='/net/scratch2/fjs2/transfer/160909/KvLip121_3/'
%workingDir='/fastscratch/fjs2/140625n/'
% workingDir='/Users/fred/EMWork/Hideki/140625/ExampleImages';
%workingDir='/Users/fred/EMWork/Hideki/150714/';
% workingDir='/Volumes/WD2Blue/EMWork/Hideki/150620/BoxTmp1_slot4_WalkerBx1_30sec/';


compressedDir=[workingDir 'Compressed/'];
localWorkingDir=workingDir;
logDir='Log/';

createMiFiles=0;

pars=struct;
pars.overwrite=1;
pars.useParfor=0;

pars.forceTracking=1;  % Drift Tracker
pars.writeZTiff=0;     % 2= z.tif only; 1= write both .mrc and z.tif
pars.writeGraphics=1;
pars.testSegments=[2 20; 38 inf];
pars.testSegments=[2 20; 25 inf]; % 161101 data
pars.writeStack=1;  % Dirft tracker

% for merging
%pars.defaultPixA=1.05;
%pars.defaultPixA=1.781;
pars.defaultPixA=0;
pars.searchDefoci=[1 8 ; 8 15]; % for MergeImages [1stmin 2ndMin ; 1stMax 2ndMax]
pars.doAlignment=1;  % MergeImages
pars.doFitting=0;
pars.doWriteInfos=1;
pars.weights=[1 0];  %%% single exposure
%pars.weights=1;
pars.mcDS=1;
%pars.mergeMode=3;  %%% no phase flip!
pars.mergeMode=1;   %%% normal
pars.mapMode='Kv'; 

pars.UsePWFilter=doPrelimInverseFilter;
pars.doPreSubtraction=1;  % rsRefineVesicleFits: pre-subtract old vesicle fit.
pars.loadFilenames=1; % pick up allNames.mat in base directory
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
    allNames=f2FindInfoFiles;
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

jobNames=allNames(ourBlockStart:ourBlockEnd);
numJobNames=numel(jobNames);
if numJobNames<1  % nothing to do
    error(['No mi files found: ' pwd]);
end;

if serialMode
    nNames=1;
else
    nNames=numJobNames;
end;
iName=1:nNames;
while iName(end)<=numJobNames
    disp(['Working on images ' num2str(iName(1)) ' to ' num2str(iName(end)) ...
        ' of ' num2str(numJobNames)]);

    ourNames=jobNames(iName);
    if serialMode
        disp(ourNames{1});
        mi=ReadMiFile(ourNames{1});
        [logSequence,dates]=miDecodeLog(mi);
    else
        logSequence=true(1,10);
    end;
   
    % find jump (sequence 1)
    if doFindJump && ~f2Mode
        k2FindDefocusJump(ourNames,pars);
    end;
   
    % drift tracker (sequence 2)
    if doTrack
        if f2Mode
            f2DriftTracker(ourNames,pars);
        else
            k2DriftTracker(ourNames,pars);
        end;
    end;
    % merge images (sequence 3)
    if doMerge
        if ~logSequence(3) || forceMerging || now-dates(3)>maxAge
            MergeImages(ourNames,pars);
        else
            disp('Merging skipped.');
        end;
    end;

    doCompressMovies=doCompressMovies && f2Mode;  % can't do this with k2 movies.
    if doCompressMovies || doCompressMicrographs
        f2CompressMovies(ourNames,compressedDir,doCompressMovies,doCompressMicrographs);
    end;
    if doDownsampleMergedImages
        DownsampleMergedImages(ourNames);
    end;
    % inverse filter (sequence 6)
    if doPrelimInverseFilter && logSequence(6)<=logSequence(3)
        fpars=struct;
        fpars.useUnsubImage=1;
        meInverseFilterAuto(ourNames,fpars);
    end;
    % find vesicles (sequence 4) *****************
    if doFindVesicles && (logSequence(4)<=logSequence(3) || now-dates(4)>maxAge)
            VesicleFinder(ourNames);
    elseif doMultiFindVesicles
        for i=1:numel(findVesicleAmps)
            vfpars.sav.vesicleAmps=[findVesicleAmps(i) 2e-3 0];
            cd('..');
            cd(findVesicleDirs{i});
            disp(pwd);
            VesicleFinder(ourNames,vfpars);
        end;
    end;
    % refine vesicles (sequence 5)
    if doRefineVesicles && (logSequence(5)<logSequence(4) ...
            || logSequence (5)< minRefineVesiclesSequence || now-dates(5)>maxAge)
        rpars=pars;
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
    elseif doRefineVesicles
        disp('  Refine Vesicles skipped.');
    end;
    % inverse filter (sequence 6)
    if doInverseFilter && (logSequence(6) <= logSequence(5)) ...
            || now-dates(6)>maxAge
        meInverseFilterAuto(ourNames);
    elseif doInverseFilter
        disp('  Inverse Filter skipped.');
    end;
    % picking preprocessor (sequence 8)
    if doPickingPreprocessor && (logSequence(8) <= logSequence(6))
        % no picking after latest vesicle refinement? Then run it.
       rsPickingPreprocessor4(ourNames,pars);
    end;
    
    iName=iName+nNames;
end;
