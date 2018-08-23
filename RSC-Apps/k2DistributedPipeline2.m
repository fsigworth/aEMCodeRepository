% k2DistributedPipeline.m
% Run the processing pipeline for K2 or F2 movies.
% This is assumed to run after running k2CreateInfoFiles or f2CreateInfoFiles
% on one worker so
% we only read the mi files, not create them.

f2Mode        =0;  % 0 means k2 data
serialMode    =0;

doFindJump    =0;
doTrack       =1;  % do movie alignment
doMerge       =1;
doCompressMovies      =0;  % compress movies
doCompressMicrographs =0;  % compress micrographs
doRefineVesicles=0;
doInverseFilter=0;
doPickingPreprocessor=0;

%workingDir='/ysm-gpfs/scratch60/fjs2/170609/KvLipo135_1/'
workingDir='/ysm-gpfs/scratch60/fjs2/170808/KvLipo/'
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
pars.overwrite=0;

pars.useParfor=0;
pars.forceTracking=1;  % Drift Tracker
pars.writeZTiff=1;     % 2= z.tif only; 1= write both .mrc and z.tif
pars.writeGraphics=1;
pars.testSegments=[2 20; 38 inf];
%pars.testSegments=[2 20; 45 inf]; % 161101 data

% for merging
pars.defaultPixA=1.533;
pars.searchDefoci=[0.5 5 ; 6 15]; % for MergeImages [1stmin 2ndMin ; 1stMax 2ndMax]
pars.doAlignment=1;  % MergeImages
pars.doFitting=1;
pars.doWriteInfos=1;
pars.weights=[1 1];
pars.mcDS=1;

pars.loadFilenames=1; % pick up allNames.mat in base directory
pars.cpe=0;  % 0 means no change.

pars.modelMiName='~/scratch60/170417/KvLipo134_4/sq02w11/Info/sq02_1_0001_Apr18_15.11.43mi.txt';


doSimulateBatch=0;  % for simulating batch on local machine

% Figure out who we are
host=getenv('HOSTNAME');
isCluster=strncmpi(host,'c',1);  % farnam nodes are like 'c22n09'
disp(['isCluster = ' num2str(isCluster)]);
sb=getenv('BATCH_ENVIRONMENT');
isInteractive=~( strcmp(sb,'BATCH') || doSimulateBatch );
disp(['isInteractive = ' num2str(isInteractive)]);
pars.showGraphics=isInteractive || doSimulateBatch;  % show graphics anyway if simulating.

numJobs=str2double(getenv('NUM_JOBS'));
jobIndex=str2double(getenv('JOB_ID')); % one-based index

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
if isInteractive
    nGroups=MyInput('nGroups',1);
    if nGroups>1
        iGroup=MyInput('iGroup',1);
    else
        iGroup=1;
    end;
else  % batch mode.  numNodes must be even or else we have a duplicate group.
    iGroup=jobIndex; % will be 1, 2, etc.
    nGroups=numJobs;
end;

if isCluster
    cd(workingDir);
else
    cd(localWorkingDir)
end;

% Set up the log file
CheckAndMakeDir('Log',1);

% Create a log file.
logName=[logDir sprintf('%02d',iGroup) 'log.txt'];
logFile=fopen(logName,'a');
logs=struct;
logs.handles=[1 logFile];
logs.idString='';
pars.logs=logs;
mdisp(pars.logs,' ');
mdisp(pars.logs,'======================================================================-');
mprintf(pars.logs.handles,'%s%s%d\n',datestr(now),' Startup: group ',iGroup);
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

ngx=2;

% Figure out what we will do
blockSize=(nNames-1)/nGroups;
% group 1 has an extra ngx entries.
ourBlockStart=round((iGroup-1)*blockSize)+1+ngx*(iGroup>1);
ourBlockEnd=min(nNames,round(iGroup*blockSize)+ngx);
mprintf(pars.logs.handles,'files %d to %d\n',ourBlockStart,ourBlockEnd);

jobNames=allNames(ourBlockStart:ourBlockEnd);
nJob=numel(jobNames);
if nJob<1  % nothing to do
    error(['No mi files found: ' pwd]);
end;

if serialMode
    nNames=1;
else
    nNames=nJob;
end;
iName=1:nNames;
while iName(end)<=nJob
    ourNames=jobNames(iName);
    if doFindJump && ~f2Mode
        k2FindDefocusJump(ourNames,pars);
    end;
    if doTrack
        if f2Mode
            f2DriftTracker(ourNames,pars);
        else
            k2DriftTracker(ourNames,pars);
        end;
    end;
    if doMerge
        MergeImages(ourNames,pars);
    end;
    doCompressMovies=doCompressMovies && f2Mode;  % can't do this with k2 movies.
    if doCompressMovies || doCompressMicrographs
        f2CompressMovies(ourNames,compressedDir,doCompressMovies,doCompressMicrographs);
    end;
    
    if doRefineVesicles
        rsRefineVesicleFits(ourNames,pars);
    end;
    if doInverseFilter
        meInverseFilterAuto(ourNames);
    end;
    if doPickingPreprocessor
       rsPickingPreprocessor4(ourNames);
    end;
    iName=iName+nNames;
end;
