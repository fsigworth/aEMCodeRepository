% k2DistributedVesFind.m
% Run the processing pipeline for Vesicle_finding_GUI and rsRefineVesicleFits
% (This is a simplified version of k2DistributedPipeline)

progName='k2DistributedVesFind.m';

doFindVesicles = 0;
doRefineVesicles      =0;
forceRefineVesicles   =0;
refineVesicleAmpsOnly=0;
doPickingPreprocessing =1;

maxSkipAge=1;  % if the corresponding log entry has a date stamp < maxAge
% days before the present we dont process it.

% workingDir='/gpfs/ysm/scratch60/sigworth/fjs2/20201203/'
workingDir='/gpfs/ysm/scratch60/sigworth/hs468/DataFromRIKEN/201228/025035_1_1/';
infoDir='Info_6/';
logDir='~/Logs_6/';

pars=struct;

pars.doPreSubtraction=1;  % rsRefineVesicleFits: pre-subtract old vesicle fit.
pars.rTerms=[100 150 200 300  inf];

pars.dsSmall=4; % downsampling of 'small' merged image

pars.loadFilenames=1; % pick up allNames.mat in base directory
pars.cpe=0;  % 0 means no change.
    pars.modelMiName='~/data/MembraneRef/160909_sq02_1_01mi.txt';


doSimulateBatch=0;  % for simulating batch on local machine
CheckAndMakeDir([workingDir 'Merged/'],1);


% Figure out who we are
host=getenv('HOSTNAME')
%isCluster=strncmpi(host,'farnam',6);  % farnam node
isCluster=any(host(1)==['c' 'p'])
isCluster=1;  

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

% Set up the log directory
CheckAndMakeDir(logDir,1);

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
mdisp(pars.logs.handles,progName);
mdisp(pars.logs,pwd);

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

jobNames=allNames(ourBlockStart:ourBlockEnd);
numJobNames=numel(jobNames);
if numJobNames<1  % nothing to do
    error(['No mi files found: ' pwd]);
end;

if doFindVesicles
    disp('------Finding Vesicles------');
    Vesicle_finding_GUI(jobNames);
end;

if doRefineVesicles
    disp(' ');
    disp('------Refining vesicles------');
    disp(' ');
    iName=1;
    while iName<=numJobNames
        %     disp(['Working on image ' num2str(iName) ' of ' num2str(numJobNames)]);
        mprintf(pars.logs.handles,'Working on image %4d of %4d\n',iName,numJobNames);
        ourName=jobNames(iName);
        disp(ourName{1});
        mi=ReadMiFile(ourName{1});
        [logSequence,dates]=miDecodeLog(mi);
        
        % refine vesicles (sequence 5)
        if forceRefineVesicles || (now-dates(5)>maxSkipAge)
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
            
            rsRefineVesicleFits(ourName,rpars);
        elseif ~forceRefineVesicles
            mprintf(pars.logs.handles,' --skipped (recently done)\n');
        end;
        
        iName=iName+1;
    end;
end;

if doPickingPreprocessing
    disp(' ');
    disp('----------Picking Preprocessing--------');
    rsPickingPreprocessor4(jobNames);
end;

mdisp(pars.logs.handles,['Finished. ' datestr(now)]);
fclose(logFile);

