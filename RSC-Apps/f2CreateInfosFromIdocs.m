function miNames=f2CreateInfosFromIdocs(ciPars)
% function miNames=f2CreateInfosFromIdocs(ciPars)
% nargin=0;
% Find .idoc files and create the 
% corresponding mi (micrograph info) files for a set of .tif summed image
% files. Parameters are read from the optional ciPars structure.  The
% output is a cell array of the mi file names.
%
% Directory structure
%   The existing .tif images and .idoc files
%   are assumed to reside in a directory like this:
% experiment/Exposures/sq01/sq01_1.idoc and
% experiment/Exposures/sq01/sq01_10000.tif, etc.
% In batch mode (either matlab was started with -nojvm, or
% ciPars.batchMode=1), we assume that we're in the experiment/ directory and
% this program looks for a directory whose name starts with 'Expos'.  In
% interactive mode, a file selector is put up and you would select the
% 'Exposures' directory.
% Either way, the subdirectories inside 'Exposures' are searched for
% files with the extension .idoc.  These are read to find the image files,
% which are assumed to be in the same directory.
% The program creates a directory Info/ at the same level as the
% Exposures directory.  Into this we write an mi file corresponding to
% each movie.  The name is constructed from the name of the first .tif file
% of the image set.  Thus in the example above we will have
%   Info/sq01_10000.txt, etc.
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
lp.batchMode=0;
lp.overwrite=0;
lp.cameraIndex=6;  % 6 means Falcon2 camera.
% Parameters that are *not* read from the idoc files:
lp.cpe=60;        % counts/electron
lp.kV=200;
lp.imageSize=[4096 4096];

idocExtension={'.idoc'};
dirInfo='Info/';

% Overwrite these with values from the ciPars struct if present.
if nargin<1
    ciPars=struct;  % default is an empty struct
end;
lp=SetOptionValues(lp,ciPars);

lp.batchMode=lp.batchMode || ~usejava('desktop');  % Matlab's gui is not active, must use batch mode


% Have the user select the movie folder, or else find it.
if ~lp.batchMode  % Matlab's GUI is active
    disp('Getting the framesPath');
    framesPath=uigetdir('.','Select a directory containing directories of image files');
    if isnumeric(framesPath)
        return
    end;
    [rootPath,expoDir]=ParsePath(framesPath);
    cd(rootPath);
else
    d=dir;
    expoDir=[];
    for i=1:numel(d)
        if strncmpi(d(i).name,'Expos',5) && d(i).isdir
            expoDir=AddSlash(d(i).name);
            break
        end;
    end;
    if numel(expoDir)<1
        error('Couldn''t find an Exposures directory');
    end;
    rootPath=AddSlash(pwd);
end;

d=dir(expoDir);  % scan for subdirectories
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
if ~exist(dirInfo,'dir')
    disp(['Creating ' dirInfo])
    mkdir(dirInfo);
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
mi0.infoPath=dirInfo;
mi0.frameDose=[];

% Loop over subdirectories
% We pick up a list mNames of idoc names
for paIndex=1:numel(pathName)
    squareDir=pathName{paIndex};
    disp(' ');
    disp(['Directory: ' expoDir squareDir]);
    d=dir([expoDir squareDir]);
    numFiles=0;
    sNames={};  % names in subdirectory
    for i=1:numel(d)  % scan the subdirectory and get idoc file names
        [~,nm,ex]=fileparts(d(i).name);
        if any(strcmpi(ex,idocExtension)) && nm(1)~='.' % we take .idoc files
            numFiles=numFiles+1;
            sNames{numFiles}=d(i).name;
        end;
    end;
    
    disp([num2str(numFiles) ' idoc files to process']);
    names=cell(0,1);
    pos=zeros(0,2);
    defoci=zeros(0,1);
    pixA=zeros(0,1);
    nImgs=0;
    doses=zeros(0,1);
    
    for findex=1:numFiles
        fname=sNames{findex};
        disp(['reading ' fname]);
        idoc=fopen([expoDir squareDir fname]);
        
        %%
        while ~feof(idoc) % scan the entire file
            s=fgetl(idoc);
            while ~feof(idoc) && numel(s)==0  % skip blanks
                s=fgetl(idoc);
            end;
            % Find an imge block
            ok=0;
            while ~feof(idoc) && ~ok
                if numel(s)>0
                    c=textscan(s,'[Image = %s','emptyvalue',0);
                    tok=char(c{1});
                    ok=numel(tok)>0;
                end;
                s=fgetl(idoc);
            end;
            if ok  % we got a name
                name=char(c{1});
                name(end)=[];  % delete the final ]
                nImgs=nImgs+1;
                names{nImgs}=name;                
                while ~feof(idoc) && numel(s)>0  % keep going till we hit a blank line.
                    c=textscan(s,'%s = %f %f %f','emptyvalue',0);
                    tok=char(c{1});
                    switch tok
                        case 'StagePosition'
                            pos(nImgs,:)=[c{2} c{3}];
                        case 'ImageShift'
                            pos(nImgs,:)=pos(nImgs,:)+[c{2} c{3}];
                        case 'TargetDefocus'
                            defoci(nImgs)=c{2};
                        case 'PixelSpacing'
                            pixA(nImgs)=c{2};
                        case 'ExposureDose'
                            doses(nImgs)=c{2};
                    end;
                    s=fgetl(idoc);
                end;  % encountered a blank line
            end;
        end;
        fclose(idoc);
    end; % for findex
end; % for paIndex
%%
% Now create mi files
totalNFiles=0;
miNames=cell(0,1);
% Construct the output name
i=1;
while i<=nImgs
    mi=mi0;
    [~,baseName,~]=fileparts(names{i});
    mi.baseFilename=baseName;
    mi.imagePath=[AddSlash(expoDir) AddSlash(squareDir)];
    mi.pixA=pixA(i);
    mi.imageFilenames{1}=names{i};
    mi.doses=doses(i);
    mi.identifier=rand;
    oldPos=pos(i,:);
    j=1;
    i=i+1;
    while i<=nImgs && all(abs(pos(i,:)-oldPos)<.1)  % add more images
        j=j+1;
        mi.imageFilenames{j}=names{i};
        mi.doses(j)=doses(i);
        i=i+1;
    end;
    
    %         Write the mi file, if desired.
    miBasename=[dirInfo mi.baseFilename 'mi'];
    miTextName=[miBasename '.txt'];
    savedName=miTextName;  % default
    str='';
    miExists=exist(miTextName,'file');  % don't already have a mi.txt file
    if ~miExists
        miMatName=[miBasename, '.mat'];
        if exist(miMatName,'file')
            mi=ReadMiFile(miMatName);
            WriteMiFile(mi,miTextName);  % convert an existing .mat to .txt
            str='Converted to text';
            miExists=true;
        end;
    else
        str='Already exists';
    end;
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
end;

