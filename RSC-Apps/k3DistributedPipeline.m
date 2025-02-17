% k3DistributedPipeline.m
% Run the processing pipeline for K2 or K3 movies.
% Simplified version does vesicle finding, refinement and picking
% preprocessor.
doFindVesicles        =1;
doPrelimInverseFilter =0;
doRefineVesicles      =0;
doInverseFilter       =0;
doPickingPreprocessor =0;

jobStart=1001;
jobEnd=inf;

% Directories must be set up
workingDir='/gpfs/ysm/scratch60/sigworth/fjs2/20211122/';
%workingDir='/gpfs/ysm/scratch60/sigworth/fjs2/20220125/';
cd(workingDir);

%workingDir='/Users/fred/EMWork/Yangyu/20211122_sel/';
infoDir='Info_C24-4/';
% infoDir='Info_C15-2/';
%infoDir='Info/';
logDir='Log/';

pars=struct;
pars.loadFilenames=0; % pick up allNames.mat in base directory
pars.filenameFile='allNames.mat'; % in the working directory

% for picking preprocessor
pars.mapMode='Kv';

% ---for rsRefineVesicleFits----
rpars=struct;
%  rpars.doPreSubtraction=1;  % rsRefineVesicleFits: pre-subtract old vesicle fit.
%  rpars.rTerms=[100 150 200 300  inf];
%  rpars.rTerms=[90 100 120 150 200 250 300 inf];
rpars.peakPositionA=[];
% rpars.dsSmall=4; % downsampling of 'small' merged image
rpars.skipFitting=0;
rpars.writeSubMRC=1;
rpars.writeSmallMRC=0;
rpars.writeSmallSubMRC=1;

% Figure out who we are
host=getenv('HOSTNAME');
isCluster=any(strncmpi(host,{'c' 'p'},1));  % farnam nodes are like 'c22n09'
disp(['isCluster = ' num2str(isCluster)]);
sb=getenv('BATCH_ENVIRONMENT');
isInteractive=~( strcmp(sb,'BATCH'));
disp(['isInteractive = ' num2str(isInteractive)]);
pars.showGraphics=isInteractive;  % show graphics anyway if simulating.

numJobs=str2double(getenv('NUM_JOBS'));
jobIndex=str2double(getenv('JOB_ID')); % one-based index

if isInteractive
    disp('Batch simulation:');
    numJobs=MyInput('numJobs',1);
    jobIndex=MyInput('jobIndex',1);

else
    if isnan(numJobs)
        disp('undefined number of jobs');
        numJobs=1;
    end;
    if isnan(jobIndex)
        disp('undefined job index');
        jobIndex=1;
    end;
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

if pars.loadFilenames
    disp('loading filenames');
    if exist('allNames.mat','file')
        load allNames.mat
    else
    disp('allNames.mat not found. Finding the mi files');
    allNames=f2FindInfoFiles(infoDir);
    end;
else
    disp('Finding the mi files');
    allNames=f2FindInfoFiles(infoDir);
end;
nNames=numel(allNames);
disp([num2str(nNames) ' files total']);

if nNames<1
    msg=['No mi files found in ' pwd filesep infoDir];
    mdisp(pars.logs,msg);
    error(msg);
end;

% Figure out what we will do
% jobOffset tells how many files to skip at beginning or at end
jobStart=max(1,min(nNames,jobStart));
jobEnd=min(nNames, jobEnd);
blockSize=(jobEnd-jobStart+1)/numJobs;
ourBlockStart=round(jobStart+(jobIndex-1)*blockSize);
ourBlockEnd=min(jobEnd,round(jobIndex*blockSize+jobStart-1));
mprintf(pars.logs.handles,'files %d to %d\n',ourBlockStart,ourBlockEnd);

ourNames=allNames(ourBlockStart:ourBlockEnd);
numNames=numel(ourNames);
if numNames<1  % nothing to do
    error(['No mi files found: ' pwd]);
end;

    % inverse filter (sequence 6)
    if doPrelimInverseFilter
        fpars=struct;
        fpars.useUnsubImage=1;
            meInverseFilterAuto(ourNames,fpars);
    end;

    if doFindVesicles
            Vesicle_finding_GUI(ourNames);
        end;

        % refine vesicles (sequence 5)
    if doRefineVesicles
            rsRefineVesicleFitsA(ourNames,rpars);
    end;
            % inverse filter (sequence 6)
    if doInverseFilter 
            meInverseFilterAuto(ourNames);
    end;

    % picking preprocessor (sequence 8)
    if doPickingPreprocessor
            rsPickingPreprocessor4(ourNames,pars);
    end;
    