function MergeImages(miNames,mpars)
% function MergeImages(miNames,mpars)
% can also be called by
%   MergeImages(parFileName)
%   MergeImages % put up file selector and use defaults
%
% Same as the old MergeImages, but operates on already extant mi files,
% reading names from the cell array miNames or putting up a file selector.
% This is a version that works in parallel on Louise nodes.
% In batchMode both arguments are required.  For interactive mode, no
% arguments are necessary unless you want to override defaults using the
% struct mpars.
% mpars.ctfOptions is a struct that contains further options for ctf
% determination.
%
%
% Some useful mpars values (and defaults)
% mpars.overwrite (1)  Forces new merging to be performed
% mpars.initialDefoci (1.5 10] Initial values for ctf fitting
% mpars.mergeMode (1)  Changes filtering of the merged images

if nargin<2
    mpars=struct;
else
    save('MergeImagesPars.mat','mpars');
    disp('Saving MergeImagesPars');
end;
if nargin<1
    miNames={};
end;

% Set up the default parameter values.
suffixes={'m' 'msf' 'msu' 'msi'}; % Suffix for merged image name, indexed by mergeMode.
suffixes={'m' 'm' 'm' 'm'}; % Suffix for merged image name, indexed by mergeMode.
pars=struct;

% See if we are on the Louise cluster.  By default we use parfor if we are
% on the cluster.
% host=getenv('HOSTNAME');
% pars.useParfor=strncmp(host,'compute',7);
pars.useParfor=0;

mctfOptions=struct;
% if isa(miNames,'char')  % we got a string as argument
%     mdisp(pars.logs,['Reading parameters: ' miNames]);
%     mpars=ReadMiText(miNames);
    if isfield(mpars,'ctfOptions')
        mctfOptions=mpars.ctfOptions;
        mpars=rmfield(mpars,'ctfOptions');
    end;
% end;

pars.overwrite=1;
pars.batchMode=1;        % by default, we're not interactive
pars.useParfor=0;
pars.writeZTiff=0;        % 0=mrc only; 1=both; 2=z.tif only.
pars.showGraphics=1;
pars.writeGraphics=1;   % save graphics files
pars.preWhitenNoise=0;   % pre-whiten micrographs before mergint (for CCD)
pars.mcDS=2;             % Composite image downsample ratio
pars.nZeros=1;           % number of zeros limiting high-defocus contributions
pars.searchDefoci=[0.8 8 ; 5 20];
pars.weights=1;         % use the existing mi.weights unless more than one element is given.
pars.doFitting=1;     % If alignments and CTFs already determined, just merge.
pars.doAlignment=1;   % determine CTFs but not alignment
pars.doWriteInfo=1;        % Store the new mi file
pars.doMakeDirectories=1;  % if directories don't exist, create them.
pars.ignoreOldCTFs=1;    % If a CTF fit already exists, don't use it to set the defocus range.
pars.useCircMask=1;      % Apply circular antialias filter in meCombineImages
pars.makeJpegs=1;
pars.mergeMode=1;  % Normal merging with flat noise spectrum
% mergedMode=2 is phase-flipped but with 1st exposure CTF; 3 is not phase flipped.
pars.normalizeMicrographs=1;
pars.smallImageSize=960;
% mergeMode=[1 2] for old (pre 8/2016) merging process
pars.logs=1;
pars.overridePixA=0;  % set to new pixA value if desired.
pars.removeOutliers=0;

pars.mergeJPath='Merged/jpeg/';
pars.filtPath = 'Merged/filtered/';
pars.ctfJPath = 'Jpeg/';

startIndex=1;
endIndex=inf;

ctfOptions=struct;
ctfOptions.lowB=1;  % lower B factors
ctfOptions.spFromWholeMicrograph=1;

ctfOptions.maxRes=[8 15];
ctfOptions.minRes=[40 60];
ctfOptions.defocusStep=[.2 .4];
ctfOptions.B=40;
ctfOptions.alpha=.02;

pars=SetOptionValues(pars,mpars);
ctfOptions=SetOptionValues(ctfOptions,mctfOptions);

pars.suffix=suffixes{pars.mergeMode(1)};  % suffix for output image name.
useParfor=pars.useParfor

mdisp(pars.logs,'--------MergeImages--------');

% We should put up the file selector if no filenames were given
if ~isa(miNames,'cell') || numel(miNames)<1
    mdisp(pars.logs,'No mi filenames were given, opening file selector...');
    pars.batchMode=0;  % try to force the user to find filenames
end;
% But, force batch mode if we have no GUI
pars.batchMode=pars.batchMode || (~isdeployed && ~usejava('desktop'));

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
else
    rootPath=AddSlash(pwd);
end;


%% -------------- Main program ------------ %%
if pars.showGraphics
    figure(10);
    SetGrayscale;
    figure(1);
    SetGrayscale;
end;
%%

nmi=numel(miNames);
lastImage=min(nmi,endIndex);
mi=ReadMiFile(miNames{startIndex},1);
nex=numel(mi.imageFilenames);  % number of exposures

% % Set up the initial defocus values
% if nex>numel(pars.initialDefoci)  % replace with sequence (1 4 8 16....) x first.
%     pars.initialDefoci(2)=4*pars.initialDefoci(1);
%     pars.initialDefoci(3:nex)=pars.initialDefoci(2)*2.^(1:nex-2);
%     mdisp(pars.logs,['initialDefoci = ' numstr(pars.initialDefoci)]);
% end;
% 

modelSpectrum=CCDModelSpectrum2D(mi.camera);  % handle DE-12 or CCD or k2
% [nCx, nCy]=size(modelSpectrum);  % Use the non-normalized spectrum, as
% mePreWhiten normalizes it.

mi.basePath=rootPath;

% Get all the directories we'll use for writing results
writePaths={[mi.basePath mi.infoPath]
    [mi.basePath mi.procPath]
    [mi.basePath pars.mergeJPath]
    [mi.basePath pars.filtPath]
    [mi.basePath pars.ctfJPath]};

%if needed, create directories.  The imagePath must already exist.
for j=1:numel(writePaths);
    if ~DirectoryExists(writePaths{j})
        if pars.doMakeDirectories
            mkdir(writePaths{j});
        else
            warning(['the path doesn''t exist: ' writePaths{j}]);
        end;
    end;
end;

tic
mdisp(pars.logs,['Merging image sets ' num2str(startIndex) ' to ' num2str(lastImage)]);
if pars.useParfor
    parfor imIndex=startIndex:lastImage
        imIndex
        mdisp(pars.logs,['--Image index ' num2str(imIndex)]);
        mi=ReadMiFile(miNames{imIndex},1);
        mi.basePath=AddSlash(rootPath);
        %         Do the processing here
        [mi,ok]=MergeImagesFcn(mi,pars,ctfOptions,modelSpectrum);
        %
        if ok && pars.doWriteInfo
            mi.log{end+1,1}=['MergeImages ' TimeStamp];
            outName=WriteMiFile(mi,miNames{imIndex});
            mdisp(pars.logs,['Updated: ' outName]);
        elseif ok
            outName=miNames{imIndex};
            mdisp(pars.logs,['Not written: ' outName]);
        else
            outName=miNames{imIndex};
            mdisp(pars.logs,['Not processed: ' outName]);
        end;
    end;
else
    for imIndex=startIndex:lastImage
        mdisp(pars.logs,['--Image index ' num2str(imIndex)]);
        mi=ReadMiFile(miNames{imIndex},1);
        mi.basePath=AddSlash(rootPath);
        mi.weights=pars.weights;  % Overwrite mi weights
        %         Do the processing here
        [mi,ok]=MergeImagesFcn(mi,pars,ctfOptions,modelSpectrum);
        %
        if ok && pars.doWriteInfo
            mi.log{end+1,1}=['MergeImages ' TimeStamp];
            outName=WriteMiFile(mi,miNames{imIndex});
            mdisp(pars.logs,['Updated: ' outName]);
        elseif ok
            outName=miNames{imIndex};
            mdisp(pars.logs,['Not written: ' outName]);
        else
            outName=miNames{imIndex};
            mdisp(pars.logs,['Not processed: ' outName]);
        end;
    end;
 end;
toc

nSets=lastImage-startIndex+1;
if nSets==1
    mdisp(pars.logs,'1 image set processed.');
else
    mdisp(pars.logs,[num2str(nSets) ' image sets processed.']);
end;
mdisp(pars.logs,' ');
