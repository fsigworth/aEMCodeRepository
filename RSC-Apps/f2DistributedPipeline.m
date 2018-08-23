% f2DistributedPipeline.m
% This is assumed to run after running f2CreateInfoFiles on one worker.

runDriftTracker=1;
runMergeImages=1;

pars=struct;
pars.overwrite=1;
pars.useParfor=0;
pars.writeZTiff=0;

hostname=getenv('HOSTNAME');
hostname(5:end)=' '; % pad with blanks
switch hostname(1:4)
    case 'Katz'
%         workingDir='/Users/fred/EMWork/Hideki/140625/KvBetaLiposome_new_pH8_KvLipo11_slot1';
        workingDir='/Users/fred/EMWork/Hideki/150114/Testing/';
    case 'comp'
        workingDir='/fastscratch/fjs2/140625';
    case 'tita'
        workingDir='/fs/Arctica_Data/shigematsu/Arctica/SerialEM/data/150714/Tor2slot4/';
    otherwise
        workingDir='';
end;

cd(workingDir);

logDir='Log/';
doSimulateBatch=0;

% Figure out who we are

s=getenv('PBS_ENVIRONMENT');
sb=getenv('BATCH_ENVIRONMENT');
isInteractive=~(strcmp(s,'PBS_BATCH') || strcmp(sb,'BATCH') || doSimulateBatch);
pars.showGraphics=isInteractive;

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


% Set up the log directory
CheckAndMakeDir(logDir,1);

% Create a log file.
logName=[logDir sprintf('%02d',jobIndex) 'log.txt'];
logFile=fopen(logName,'a');
logs=struct;
logs.handles=[1 logFile];
logs.idString='';
pars.logs=logs;
mdisp(pars.logs,'----------------------------');
mprintf(pars.logs.handles,'%s%s%d\n',datestr(now),' Startup: node ',jobIndex);

allNames=f2FindInfoFiles;
nNames=numel(allNames);

if nNames<1
    mdisp(pars.logs,'No mi files found.');
    error('No mi files found');
end;

blockSize=ceil(nNames/numJobs);
ourBlockStart=jobIndex*blockSize+1;
ourBlockEnd=min(nNames,ourBlockStart+blockSize-1);
mprintf(pars.logs.handles,'files %d to %d\n',ourBlockStart,ourBlockEnd);

ourNames=allNames(ourBlockStart:ourBlockEnd);
if numel(ourNames)<1  % nothing to do
    return
end;

if runDriftTracker
    f2DriftTracker(ourNames,pars);
end;
if runMergeImages
    MergeImages(ourNames,pars);
end;
