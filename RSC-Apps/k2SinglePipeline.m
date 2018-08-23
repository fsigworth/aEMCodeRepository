% k2SinglePipeline.m

% Assume we've cd'd to the experiment folder.

pars=struct;
pars.getAllNames=1;
pars.kV=200;
pars.batchMode=1;
pars.defaultPixA=1.247;
pars.nFrameSets=2;
pars.readOnly=0;
pars.useParfor=0;


doSimulateBatch=0;

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
    cd('/fastscratch/fjs2/151117/KvLipo50slot3');
else
    cd('/Users/fred/EMWork/Hideki/151117/KvLipo50slot3');
end;

phase=jobIndex;
stride=numJobs;



while 1
    
    names=k2CreateInfoFiles(pars);
    pause(1);  % Wait for incomplete files to be written
    
    pars.readOnly=0;
    nim=numel(names);
    %%
    iStart=nim;  % work backwards
    iStart=iStart-stride+phase+1;
    for i=iStart:-stride:1
        % for i=iStart:nim
        disp(['processing ' names{i}]);
        k2FindDefocusJump(names(i),pars);
        k2DriftTracker(names(i),pars);
        MergeImages(names(i),pars);
    end;
    
end;