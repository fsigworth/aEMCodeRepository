function miNames=k2FindDefocusJump(miNames,fdPars)

% Attempt to find the frame where the defocus jump occurs in a K2 movie.
% Some useful parameters (and their default values) are
% fdPars.batchMode (true)
% fdPars.testSegments ([2 17; 25 inf])
% fdPars.overwrite (false)
% fdPars.initialDefoci ([1.5 10])
%
% Sets the following fields in the mi file:
% frameDose
% frameSets
% In case the user used the file selector, the new miNames is returned.


% Set overall defaults, which can be overriden by corresponding fdPars fields.
pars=struct;

% See if we are on the Louise cluster, and by default use parfor
host=getenv('HOSTNAME');
% Here are the defaults for all the parameters.
% the most important are
% testSegments and initialDefoci
pars.useParfor=strncmp(host,'compute',7);  % by default, use parfor on the cluster
pars.batchMode=true;
pars.testSegments=[2 35; 45 inf];  % note: testSegments(1,1) sets first frame to keep!
pars.searchDefoci=[.8 5 ; 6 12];
pars.overwrite=1;
pars.writeGraphics=1;
pars.showGraphics=1;
pars.ds=8;
pars.ds1=8;
pars.maxRes=[10 20];
pars.minRes=[40 80];
pars.kV=200;
pars.dirJpeg='Jpeg/';
pars.logs=1;  % default log

if nargin<1
    miNames={};
end;
if nargin<2
    fdPars=struct;
end;
pars=SetOptionValues(pars,fdPars);

ctfOptions.lowB=1;  % lower B factors
ctfOptions.spFromWholeMicrograph=1;
ctfOptions.defocusStep=[.3 .6];
ctfOptions.B=80;
ctfOptions.alpha=.02;

mdisp(pars.logs,'---k2FindDefocusJump---');
% We should put up the file selector if no filenames were given
if numel(miNames)<1
    mdisp(pars.logs,'No mi filenames were given...');
    pars.batchMode=0;  % try to force the user to find filenames
end;
% But, force batch mode if we have no GUI
pars.batchMode=pars.batchMode || ~usejava('desktop');

if  ~pars.batchMode
    mdisp(pars.logs,'Asking for info file names.');
    % Have the user select some info files, and copy these into miNames
    [names, pathName]=uigetfile({'*mi.txt';'*mi.mat'},'Select mi files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(names)
        names={names};
    end;
    [rootPath,infoPath]=ParsePath(pathName);
    cd(rootPath);
    miNames=cell(0);
    for i=1:numel(names)
        miNames{i}=[infoPath names{i}];
    end;
end;

nFiles=numel(miNames);
nTestSegs=size(pars.testSegments,1);

    if ~exist(pars.dirJpeg,'dir')
            mkdir(pars.dirJpeg);
    end;


if nFiles>1
    mdisp(pars.logs,['Working on ' num2str(nFiles) ' files.']);
end;
nSteps=size(pars.testSegments,1)-1;
transFrames=zeros(nFiles,nSteps);
%%
if pars.useParfor
    parfor ind=1:nFiles
        infoName=miNames{ind};
        miStart=ReadMiFile(infoName,1);
        disp(infoName);
        nFrameSets=size(miStart.frameSets,1);
        if pars.overwrite || nFrameSets<nTestSegs
            mi=k2FindJumpFcn(miStart,pars,ctfOptions);
            mdisp(pars.logs,'frameSets:');
            mdisp(pars.logs,mi.frameSets)
            mi.log{end+1}=['k2FindDefocusJump ' TimeStamp];
            fullName=WriteMiFile(mi,infoName);
            %%
            mdisp(pars.logs,['Updated: ' fullName]);
            if size(mi.frameSets,1)>1
                transFrames(ind,:)=mi.frameSets(1:end-1,2)'+1;
            end;
            mdisp(pars.logs,' ');
        else
            mdisp(pars.logs,['Skipped:  ' infoName]);
            mdisp(pars.logs,' ');
        end;
    end;
    
else
    
    for ind=1:nFiles
        infoName=miNames{ind};
        miStart=ReadMiFile(infoName,1);
        mdisp(pars.logs,infoName);
        % miStart.frameSets
        nFrameSets=size(miStart.frameSets,1);
        if pars.overwrite || nFrameSets<nTestSegs
            movieName=[miStart.moviePath miStart.movieFilename];
            if exist(movieName,'file')
                mi=k2FindJumpFcn(miStart,pars,ctfOptions);
                mdisp(pars.logs,'frameSets:');
                mdisp(pars.logs,mi.frameSets);
                % mi.frameSets
                mi.log{end+1,1}=['k2FindDefocusJump ' TimeStamp];
            else
                mdisp(pars.logs,['---Movie not found for ' infoName]);
                system(['rm ' infoName]);
                infoName=[infoName(1:end-4) '_x' infoName(end-3:end)];
                mi=miStart;
            end;
            fullName=WriteMiFile(mi,infoName);
            %%
            mdisp(pars.logs,['Updated: ' fullName]);
            if size(mi.frameSets,1)>1
                transFrames(ind,:)=mi.frameSets(1:end-1,2)'+1;
            end;
            mdisp(pars.logs,' ');
        else            
            mdisp(pars.logs,['Skipped:  ' infoName]);
            mdisp(pars.logs,' ');
        end;
    end;
end;
%
% Show the final results
if pars.showGraphics
    nr=3;
    nc=4;
    subplot(nr,nc,nr*nc-1);
    hist(transFrames);
    title('Transition frames');
end;
mdisp(pars.logs,' ');
mdisp(pars.logs,'Summary of transition frames:');
for i=1:nFiles
    mdisp(pars.logs,[miNames{i} '   ' num2str(transFrames(i,:)')]);
end;
