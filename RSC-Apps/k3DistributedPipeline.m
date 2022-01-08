% k3DistributedPipeline.m
% Run the processing pipeline for K2 or K3 movies.
% Simplified version does vesicle finding, refinement and picking
% preprocessor.

doFindVesicles        =1;
doPrelimInverseFilter =0;
doRefineVesicles      =0;
doInverseFilter       =0;
doPickingPreprocessor =0;

% Directories must be set up
workingDir='/gpfs/ysm/scratch60/sigworth/fjs2/20211122/';

workingDir='/Users/fred/EMWork/Yangyu/20211122_sel/';
infoDir='Info_all_temp/';
logDir='Log/';

pars=struct;
pars.loadFilenames=0; % pick up allNames.mat in base directory
pars.filenameFile='allNames.mat'; % in the working directory

% for picking preprocessor
pars.mapMode='Kv';

% ---for rsRefineVesicleFits----
pars.doPreSubtraction=1;  % rsRefineVesicleFits: pre-subtract old vesicle fit.
% pars.rTerms=[100 150 200 300  inf];
pars.rTerms=[90 100 120 150 200 250 300 inf];
% pars.dsSmall=4; % downsampling of 'small' merged image
% pars.overwrite=1;
pars.modelMiName='~/data/MembraneRef/160909_sq02_1_01mi.txt';

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

ngx=0; % extra names to worker 1

% Figure out what we will do
blockSize=(nNames)/numJobs;
% group 1 has an extra ngx entries.
ourBlockStart=round((jobIndex-1)*blockSize)+1+ngx*(jobIndex>1);
ourBlockEnd=min(nNames,round(jobIndex*blockSize)+ngx);
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
        rpars=struct;
%                 rpars.fitModes={'LinOnly'};
%                 rpars.fractionStartingTerms=1; % total terms to use in each round
%                 rpars.fractionAmpTerms=1;
%                 % Extra peaks in the scattering profile
%                 rpars.peakPositionA=[-37 0 37];  % empirical default.  Works a bit better than [37 37]
%                 rpars.targetPixA=10;  % downsampled image resolution for radius fitting
%                 
%                 rpars.xPeakSigmaA={5 5}; % width of extra Gaussian peaks, in angstrom
%                 %     The following must have at least as many elements as dpars.fitModes!            
            rsRefineVesicleFits(ourNames,rpars);
    end;
            % inverse filter (sequence 6)
    if doInverseFilter 
            meInverseFilterAuto(ourNames);
    end;

    % picking preprocessor (sequence 8)
    if doPickingPreprocessor
            rsPickingPreprocessor4(ourNames,pars);
    end;
    