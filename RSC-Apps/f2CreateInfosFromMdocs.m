% function miNames=f2CreateInfosFromMdocs(ciPars)
% function miNames=f2CreateInfosFromMdocs(ciPars)
% derived from f2CreateInfosFromIdocs.  ---not done yest
nargin=0;
% Find .mdoc files and create the
% corresponding mi (micrograph info) files for a set of .mrc extracted image
% files. Parameters are read from the optional ciPars structure.  The
% output is a cell array of the mi file names.
%
% Directory structure
%   The existing .mrc stack and .mdoc files
%   are assumed to reside in a directory like this:
% Exposure/sq01/sq01_1.mrc.mdoc and
% Exposure/sq01/sq01_1.mrc, etc.
% In batch mode (either matlab was started with -nojvm, or
% ciPars.batchMode=1), we assume that we're in the Exposure/ directory and
% this program looks for a directory whose name starts with 'Ser'.  In
% interactive mode, a file selector is put up and you would select the
% 'Exposure' directory.
% Either way, the subdirectories inside 'Exposure' are searched for
% files with the extension .mdoc.  These are read to find the image stacks,
% which are assumed to be in the same directory.
% The program creates a directory Info/ at the same level as the
% SerialEM directory.  Into this we write an mi file corresponding to
% each movie.  The name is constructed from the name of the .mrc file
% of the image set.  Thus in the example above we will have
%   Info/sq01_1mi.txt, etc.
% The program also creates a directory Micrograph/ at the same level as Info/
% and unpacks the .mrc stack into files in this directory.
%
% Into each mi file are stored the relevant filenames and paths, and also
% values are set for (defaults are in parentheses)
% gp.camera     (6)  % Arctica Falcon II
% gp.cpe        (50)
% gp.kV         (200)
% gp.defaultPixA       (1.247)
%   Note that the defaultPixA is used only in the case (generic .tif) that
%   no pixA value is stored.
% The names of the mi files are stored in the cell array gp.miNames.
% Typical values are
% gp.miNames{1}='Info/sq05_001_Jun10_00.01.8mi.txt'
% gp.miNames{2}='Info/sq05_002_Jun10_00.01.38mi.txt'


disp('f2CreateInfoFiles');

% Set defaults for local parameters
lp.nFocal=1;
lp.mrcMode=1;  % write as int16.  2: float
lp.batchMode=0;
lp.overwrite=1;
lp.overwriteImages=0;
lp.cameraIndex=6;  % 6 means Falcon2 camera.
% Parameters that are *not* read from the idoc files:
lp.cpe=60;        % counts/electron
lp.kV=200;
lp.imageSize=[4096 4096];

mdocExtension={'.mdoc'};
dirExpoName='Exposure/';
dirExpoSearch=dirExpoName(1:4);
dirInfoName='Info/';
dirMicroName='Micrograph/';
dirJpegName='Jpeg/';
% Overwrite these with values from the ciPars struct if present.
if nargin<1
    ciPars=struct;  % default is an empty struct
end;
lp=SetOptionValues(lp,ciPars);

lp.batchMode=lp.batchMode || ~usejava('desktop');  % Matlab's gui is not active, must use batch mode


% Have the user select the experiment folder, or else find it.
if ~lp.batchMode  % Matlab's GUI is active
    disp(['Getting the ' dirExpoName ' Path']);
    expPath=uigetdir('.',['Select the ' dirExpoName ' folder containing folders of image files']);
    if isnumeric(expPath)
        return
    end;
    [rootPath,expDir]=ParsePath(expPath);
    cd(rootPath);
else
    d=dir;
    expDir=[];
    for i=1:numel(d)
        if strncmpi(d(i).name,dirExpoSearch,numel(dirExpoSearch)) && d(i).isdir
            expDir=AddSlash(d(i).name);
            break
        end;
    end;
    if numel(expDir)<1
        error(['Couldn''t find ' dirExpoName ' directory']);
    end;
    rootPath=AddSlash(pwd);
end;

d=dir(expDir);  % scan for directories inside "Exposure/"
nd=0;
pathName=[];
for i=1:numel(d)
    if d(i).isdir && d(i).name(1)~='.'
        nd=nd+1;
        pathName{nd}=AddSlash(d(i).name);
    end;
end;

if nd>0  % subdirectories found
    disp([num2str(nd) ' directories found']);
else
    error('No enclosed directories found');
end;

%
if ~exist(dirInfoName,'dir')
    disp(['Creating ' dirInfoName])
    mkdir(dirInfoName);
end;
if ~exist(dirMicroName,'dir')
    disp(['Creating ' dirMicroName])
    mkdir(dirMicroName);
end;
if ~exist(dirJpegName,'dir')
    disp(['Creating ' dirJpegName])
    mkdir(dirJpegName);
end;

%% ------actual work is done here-------

mi0=meCreateMicrographInfoStruct14;  % Get the mi file template
% Put in the invariant parameters
mi0.kV=lp.kV;
mi0.imageSize=lp.imageSize;
mi0.camera=lp.cameraIndex;
mi0.cpe=lp.cpe;
mi0.originalBasePath=AddSlash(rootPath);
mi0.basePath=mi0.originalBasePath;
mi0.infoPath=dirInfoName;
mi0.frameDose=[];
mi0.weights=ones(1,lp.nFocal);
totalNFiles=0;
miNames={};

% Loop over subdirectories
for paIndex=1:numel(pathName)  % loop over directories inside Exposure/
    squareDir=pathName{paIndex};
    disp(' ');
    disp(['Directory: ' expDir squareDir]);
    d=dir([expDir squareDir]);
    numFiles=0;
    sNames={};  % names in subdirectory
    for i=1:numel(d)  % scan the subdirectory and get mdoc file names
        [~,nm,ex]=fileparts(d(i).name);
        if any(strcmpi(ex,mdocExtension)) && nm(1)~='.' % we take .mdoc files
            numFiles=numFiles+1;
            sNames{numFiles}=d(i).name;
        end;
    end;

    %%
    disp([num2str(numFiles) ' mdoc files to process']);
    
    for findex=1:numFiles
        fname=sNames{findex};
        disp(['reading ' fname]);
        mdoc=fopen([expDir squareDir fname]);
        
        %%
        %  -------Find the global entries-------
        %  first, skip any blank lines
        s=fgetl(mdoc);
        while ~feof(mdoc) && numel(s)==0  % skip blanks
            s=fgetl(mdoc);
        end;
        %     Read the global parameters
        while ~feof(mdoc) && numel(s)>0  % keep going till we hit a blank line.
            c=textscan(s,'%s = %s %s','emptyvalue',0);
            tok=char(c{1});
            switch tok
                case 'PixelSpacing'
                    mi0.pixA=str2double(c{2});
                case 'ImageFile'
                    stackFilename=char(c{2});
                case 'ImageSize'
                    mi0.imageSize=str2double([c{2} c{3}]);
            end;
            s=fgetl(mdoc);
            disp(s)
        end;
        
        %         Read in the image stack
        disp(['Reading ' stackFilename]);
        stack=ReadMRC([expDir squareDir stackFilename]);
        
        pos=zeros(1,2);
        defoci=0;
        doses=0;
        
        while ~feof(mdoc)  % loop over the rest of the file
            
            %         Get the [ZValue = n] line to get the image index.
            s=fgetl(mdoc);
            while ~feof(mdoc) && ~strncmp(s,'[ZV',3)
                s=fgetl(mdoc);
            end;
            c=textscan(s,'%s = %f','emptyvalue',0);
            nImgs=c{2}+1;
            %             Scan the rest of the block
            s=fgetl(mdoc);
            while ~feof(mdoc) && numel(s)>0  % keep going till we hit a blank line.
                c=textscan(s,'%s = %f %f %f','emptyvalue',0);
                tok=char(c{1});
                switch tok
                    case 'StagePosition'
                        pos(nImgs,:)=[c{2} c{3}];
                    case 'ImageShift'
                        pos(nImgs,:)=pos(nImgs,:)+[c{2} c{3}];
                    case 'Defocus'
                        defoci(nImgs)=max(0.2,-c{2});
                    case 'ExposureDose'
                        doses(nImgs)=c{2};
                end;
                s=fgetl(mdoc);
            end;  % encountered a blank line
        end;
        %         We now have the info from an entire mdoc file.  Let's make the mi
        %         files.
        for i=1:lp.nFocal:nImgs
            mi=mi0;
            [~,baseName,~]=fileparts(stackFilename);
            mi.baseFilename=sprintf('%s_%03d',baseName,i);
            mi.imagePath=[AddSlash(dirMicroName)];
            mi.identifier=rand;
            mi.doses=doses(i:i+lp.nFocal-1);
            for j=1:lp.nFocal
                fname=sprintf('%s_%03d',baseName,i+j-1);
                mfName=[dirMicroName fname '.mrc'];
                if ~exist(mfName,'file') || lp.overwriteImages
                    WriteMRC(stack(:,:,i+j-1),mi.pixA,mfName,lp.mrcMode);
                    WriteJpeg(stack(:,:,i+j-1),[dirJpegName fname '.jpg'],1e-3);
                end;
                mi.imageFilenames{j}=[fname '.mrc'];
                mi.ctf(j).defocus=defoci(i+j-1);
            end;
            
            %         Write the mi file, if desired.
            miBasename=[dirInfoName mi.baseFilename 'mi'];
            miTextName=[miBasename '.txt'];
            savedName=miTextName;  % default
            str='';
            miExists=exist(miTextName,'file');  % don't already have a mi.txt file
            if miExists
                if lp.overwrite
                    savedName=WriteMiFile(mi,1);
                    str=[str '; overwritten'];
                end;
            else
                savedName=WriteMiFile(mi,1);
                str='Written';
            end;
            totalNFiles=totalNFiles+1;
            miNames{totalNFiles,1}=savedName;
            disp([str ': ' savedName]);
        end;  % for i
        fclose(mdoc);
    end; % for findex
end; % for paIndex
%%
