function miNames=k2DriftTracker(miNames,dtPars)
% function k2DriftTracker(miNames,dtPars)
% Like the old k2DriftTracker, but reads info files from the cell array
% of names miNames.
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

disp('k2DriftTracker');

if nargin<1
    miNames={};
end;
if nargin<2
    dtPars=struct;
end;

% Set up the main defaults, which can be overriden with dtPars
host=getenv('HOSTNAME');
% pars.useParfor=strncmp(host,'compute',7);
pars.useParfor=0;
pars.overwrite=1;
pars.batchMode=1;

% Other parameters
pars.nameAli='al';  % aligned sum suffix
pars.nameAstk='as'; % frame stack suffix
pars.doDamageComp=1;
pars.ggDamageModel=1;  % Use Grant & Grigorieff damage function
pars.nameSegment='a';
pars.forceTracking=1;     % always track shifts
pars.doRestoreImages=1;  % Compute the full-size aligned images
pars.dsw=4;            % Downsampling of working images
pars.numIters=4;
pars.fig1Size=[1000 700];  % Size of main figure, in pixels.
pars.fig1=~pars.useParfor && usejava('desktop');  % index of figure to use
pars.nFrameSets=1;
pars.searchDefoci=[0.8 8 ; 5 15];

pars.dirJpeg='Jpeg/';
pars.doSaveFigures=1;
pars.showGraphics=1;
pars.writeZTiff=0;  % write the summed stack as z.tif
pars.writeMRC=1;    % write the summed stack as .mrc
pars.writeStack=0;  % Write out the entire stack of individual frames as .mrc, no damage comp.
pars.startIndex=1;
pars.logs=1;
pars.k2Mode=1;
pars.ccCorrection='delta';  % options are 'k2' and 'delta' for removing central artifact.
pars.oldK2=0;  % 0: recent k2 data, sets rot90(3).  1: k2 data before 6/2015, sets rot90(2).
% 2: no rotation of movie.
pars.finalRot90=0; % ccw rotations applied to micrograph image, after gain reference,
% before saving it.
% Note: for OHSU k2 images, we set oldK2=1 and finalRot90=3.

pars=SetOptionValues(pars,dtPars);

% We should put up the file selector if no filenames were given
if numel(miNames)<1
    disp('No mi filenames were given...');
    pars.batchMode=0;  % try to force the user to find filenames
end;
% But, force batch mode if we have no GUI
pars.batchMode=pars.batchMode || ~usejava('desktop');

if ~pars.batchMode
    % Have the user select some info files
    [names, pathName]=uigetfile({'*mi.txt';'*mi.mat'},'Select mi files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(names)
        names={names};
    end;
    [rootPath, infoPath]=ParsePath(pathName);
    cd(rootPath);
    miNames=cell(0);
    for i=1:numel(names)
        miNames{i}=[infoPath names{i}];
    end;
end;

if ~exist(pars.dirJpeg,'dir')
    mkdir(pars.dirJpeg);
end;

numFiles=numel(miNames);
disp(['Working on ' num2str(numFiles) ' file' ess(numFiles)]);
%%
%  Set up a figure if we'll be using it.
if pars.fig1>0
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
%         disp([num2str(ind) ':  ' infoName]);
%         startMi=ReadMiFile(infoName);
%         startMi.basePath=AddSlash(pwd);
%         frameSetsPresent=numel(startMi.frameSets)>0;
%         if ~frameSetsPresent
%             startMi.frameSets=[1 inf];
%         end;
%         frameShiftsPresent=numel(startMi.frameShifts)>0;
%         if pars.overwrite || ~frameShiftsPresent
%             %%
%             mi=k2DriftTrackFcn(startMi,pars);
%             mi.basePath=AddSlash(pwd);
%             mi.log{end+1,1}=['k2DriftTracker ' TimeStamp];
%             fullName=WriteMiFile(mi,infoName);
%             disp(['Updated: ' fullName]);
%             disp(' ');
%         else
%             disp('Already aligned.');
%             disp(' ');
%         end;
%     end;
%     
% else
    
    for ind=pars.startIndex:numFiles
        infoName=miNames{ind};
        disp([num2str(ind) ':  ' infoName]);
        startMi=ReadMiFile(infoName);
        startMi.basePath=AddSlash(pwd);
        frameSetsPresent=numel(startMi.frameSets)>0;
        nFrameSets=size(startMi.frameSets,1);
        if nFrameSets>=pars.nFrameSets
            if ~frameSetsPresent
                disp('No frameSets present, aligning the entire movie.');
                startMi.frameSets=[1 inf];
            end;
            frameShiftsPresent=numel(startMi.frameShifts)>0;
            if pars.overwrite || ~frameShiftsPresent
                %             mi=k2DriftTrackFcn(startMi,pars);
                mi=f2k2DriftTrackFcn(startMi,pars);
                mi.basePath=AddSlash(pwd);
                mi.log{end+1,1}=['k2DriftTracker ' TimeStamp];
                fullName=WriteMiFile(mi,infoName);
                disp(['Updated: ' fullName]);
                disp(' ');
            else
                disp('Already aligned.');
                disp(' ');
            end;
        else
            disp(['nFrameSets is ' num2str(nFrameSets) '; skipped']);
        end;
    end;
% end;


toc

