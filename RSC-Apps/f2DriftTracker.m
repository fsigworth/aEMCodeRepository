function miNames=f2DriftTracker(names,dtPars)
% miNames={'sq05_001_Jul27_12.41.39mi.txt'};
% function f2DriftTracker(miNames,dtPars)
% Like the old k2DriftTracker, but reads info files from the cell array
% of names miNames in batch mode.
%
% Create a summed image from a K2 movie, read as an MRC file
% containing raw counts, using a gain reference if available.
% The stored aligned-summed images are padded to 3840 x 3840.
% If there are multiple segments, they are stored separately.
% In case the user used the file selector, the new miNames are returned.
%
%   Various jpegs are stored here in a Jpegs directory 2 levels up.
% experiment/Jpegs/
%
% Both the Micrograph and Jpegs directories are created if they aren't
% already present.
% Derived from DriftTrackingWholeFrame5 for the DE camera.
% F.S. 16 Sep 13
% -modified to mask the black bars in K2 micrographs when ADC channels
% fail.  31 Jan 14.


if nargin<1
    names={};
end;
if nargin<2
    dtPars=struct;
end;

% Set up the main defaults, which can be overriden with dtPars
host=getenv('HOSTNAME');
pars.useParfor=strncmp(host,'compute',7);
pars.batchMode=numel(names)>0;

pars.overwrite=0;
pars.writeZTiff=0;  % write compensated sum as z.tif
pars.writeMRC=1;    % write compensated sum as .mrc
pars.writeStack=0;  % Write out the stack of individual frames as .mrc, no damage comp.
pars.logs=1;
pars.ggDamageModel=1; % Use the Grant & Grigorieff damage model
pars.k2Mode=0;

% Other parameters
pars.nameAli='al';  % aligned sum suffix
pars.nameAstk='as'; % frame stack suffix
pars.doDamageComp=1;
pars.nameSegment='a';
pars.forceFitting=1;     % always track shifts
pars.doRestoreImages=1;  % Compute the full-size aligned images
pars.dsw=4;            % Downsampling of working images
pars.numIters=4;

pars.showGraphics=1;
pars.fig1Size=[1000 700];  % Size of main figure, in pixels.
pars.fig1=pars.showGraphics && ~pars.useParfor && usejava('desktop');  % index of figure to use

pars.dirJpeg='Jpeg/';
pars.doSaveFigures=1;
pars.startIndex=1;
pars.finalRot90=0;

pars=SetOptionValues(pars,dtPars);

mdisp(pars.logs,'f2DriftTracker');

% Force batch mode if we have no GUI
pars.batchMode=pars.batchMode || ~usejava('desktop');

if ~pars.batchMode
    % Have the user select some info files
    [names, pathName]=uigetfile('*mi.*','Select mi files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(names)
        names={names};
    end;
    [rootPath, infoPath]=ParsePath(pathName);
    cd(rootPath);
else
    infoPath='';
end;
%     Put the infoPath into the names
miNames=cell(0);

for i=1:numel(names)
    miNames{i}=[infoPath names{i}];
end;


if ~exist(pars.dirJpeg,'dir')
    mkdir(pars.dirJpeg);
end;

numFiles=numel(miNames);
mdisp(pars.logs,['Working on ' num2str(numFiles) ' files']);
%%
%  Set up a figure if we'll be using it.
if pars.showGraphics && pars.fig1>0
    figure(1);
    SetGrayscale;
    sz=get(0,'screensize');  % Root monitor size.
    pars.fig1Size=min(pars.fig1Size,sz(3:4)*0.8);  % Don't allow the figure to be too big.
    pos=sz(3:4)/2-pars.fig1Size/2;
    set(gcf,'outerposition',[pos pars.fig1Size]);  % Put the figure in the center of the root monitor.
end;

tic

% if pars.useParfor
%     parfor ind=pars.startIndex:numFiles
%         infoName=miNames{ind};
%         mdisp(pars.logs,[num2str(ind) ':  ' infoName]);
%         startMi=ReadMiFile(infoName);
%         startMi.basePath=AddSlash(pwd);
%         frameSetsPresent=numel(startMi.frameSets)>0;
%         if ~frameSetsPresent
%             startMi.frameSets=[1 inf];
%         end;
%         frameShiftsPresent=numel(startMi.frameShifts)>0;
%         if pars.overwrite || ~frameShiftsPresent
%             %%
%             [mi,ok]=f2DriftTrackFcn(startMi,pars);
%             if ok
%                 mi.basePath=AddSlash(pwd);
%                 mi.log{end+1,1}=['k2DriftTracker ' TimeStamp];
%                 fullName=WriteMiFile(mi,infoName);
%                 mdisp(pars.logs,['Updated: ' fullName]);
%                 mdisp(pars.logs,' ');
%             else
%                 mdisp(pars.logs,['Tracking error in ' infoName]);
%             end;
%         else
%             mdisp(pars.logs,'Already aligned.');
%             mdisp(pars.logs,' ');
%         end;
%     end;
%     
% else
%     
    for ind=pars.startIndex:numFiles
        infoName=miNames{ind};
        mdisp(pars.logs,[num2str(ind) ':  ' infoName]);
        startMi=ReadMiFile(infoName);
        startMi.basePath=AddSlash(pwd);
        frameSetsPresent=numel(startMi.frameSets)>0;
        if ~frameSetsPresent
            mdisp(pars.logs,'No frameSets present, aligning the entire movie.');
            startMi.frameSets=[1 inf];
        end;
        frameShiftsPresent=numel(startMi.frameShifts)>0;
        if pars.overwrite || ~frameShiftsPresent
            [mi,ok]=f2DriftTrackFcn(startMi,pars);
            if ok
                mi.basePath=AddSlash(pwd);
                mi.log{end+1,1}=['f2DriftTracker ' TimeStamp];
                fullName=WriteMiFile(mi,infoName);
                mdisp(pars.logs,['Updated: ' fullName]);
                mdisp(pars.logs,' ');
            else
                mdisp(pars.logs,['Tracking error in ' infoName]);
            end;
        else
            mdisp(pars.logs,'Already aligned.');
            mdisp(pars.logs,' ');
        end;
    end;
    
% end;

toc

