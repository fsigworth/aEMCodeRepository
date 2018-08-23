function miNames=CCDCreateInfoFiles(ciPars)
% function miNames=CCDCreateInfoFiles(ciPars)
% Find image files and create the corresponding mi (micrograph info) files,
% with parameters read from the optional ciPars structure.  The output is a
% cell array of the mi file names.
%
% Directory structure
%   The existing micrographs are assumed to reside in a directory like this:
% experiment/Micrograph/sq01_1/10000.tif.  In batch mode (either matlab
% was started with -nojvm, or
% gp.batchMode=1), we assume that we're in the experiment/ directory and
% this program looks for a directory whose name starts with 'Micrograph'.  In
% interactive mode, a file selector is put up and you would select the
% 'Micrograph' directory.
% Either way, the subdirectories inside 'Micrograph/' are searched for
% files with the extensions .mrc or .tif.  These are assumed to be raw
% images in defocus sets of 2 or 3 (set by lp.numExposures, default 2.
% Thus in the experiment/ directory, we expect to find something like
%   Micrograph/sq05/10000.tif
%   Micrograph/sq05/10001.tif
%   Micrograph/sq05/10002.tif
%               ...
%
% The program creates a directory Info/ at the same level as the
% Micrograph directory.  Into this we write an mi file corresponding to
% each set of images.  The name is constructed from the subdirectory name
% (if desired), the first image name, and the suffix 'mi'.
% In the example we will have
%   Info/sq05_10000mi.txt
%   Info/sq05_10002mi.txt
%               ...
%   Info/sq08_10000mi.txt
%   Info/sq08_10002mi.txt
%               ...
% Into each mi file are stored the relevant filenames and paths, and also
% values are set for (defaults are in parentheses)
% camera     (1)
% cpe        (8)
% kV         (200)
% defaultPixA (1.247)
%   Note that the defaultPixA is used only in the case (generic .tif) that
%   no pixA value is stored.
% The names of the mi files are stored in the returned cell array miNames.
% Example values are
%   miNames{1}='Info/sq05_10000mi.txt'
%   miNames{2}='Info/sq05_10002mi.txt'
% and so on.


disp('CCDCreateInfoFiles');

lp=struct;


% Set defaults for local parameters
lp.batchMode=0;
lp.overwrite=1;
lp.cameraIndex=1;  % 1 means U4000 camera.
lp.cpe=50;        % depends on the camera
lp.kV=200;
lp.defaultPixA=2.34;
lp.checkImagePixA=0;
lp.numExposures=2;
lp.includeDirectoryName=1;
lp.imageExtension={'.tif'};
dirInfo='Info/';


% %%% overriding options for the Arctica Falcon 2 %%%
lp.cameraIndex=6;
lp.cpe=90;
lp.defaultPixA=1.349;
lp.numExposures=2;


% Overwrite these with the ciPars struct values if present.
if nargin<1
    ciPars=struct;  % default is an empty struct
end;
lp=SetOptionValues(lp,ciPars);

lp.batchMode=lp.batchMode || ~usejava('desktop');  % Matlab's gui is not active, must use batch mode


% Have the user select the micrograph folder, or else find it.
if ~lp.batchMode  % Matlab's GUI is active
    disp('Getting the framesPath');
    micrPath=uigetdir('.','Select a Micrograph directory');
    if isnumeric(micrPath)
        return
    end;
    [rootPath,micrDir]=ParsePath(micrPath);
    cd(rootPath);
    disp(['Root path:  ' rootPath]);
    disp(['Image path: ' micrPath]);
else
    d=dir;
    micrDir=[];
    for i=1:numel(d)
        if strncmpi(d(i).name,'Micrograph',10) && d(i).isdir
            micrDir=AddSlash(d(i).name);
            break
        end;
    end;
    if numel(micrDir)<1
        error('Couldn''t find a micrograph directory inside the current directory');
    end;
    rootPath=AddSlash(pwd);
end;

% Find subdirectories
d=dir(micrDir);  % scan for subdirectories
nd=0;
sdirName=[];
disp('Subdirectories:');
for i=1:numel(d)
    if d(i).isdir && d(i).name(1)~='.'  % any directory except .*
        nd=nd+1;
        sdirName{nd}=d(i).name;
        disp(d(i).name);
    end;
end;
disp(' ');

if nd>0  % subdirectories found
    disp([num2str(nd) ' directories found']);
else
    error('No enclosed directories found');
    return
end;


if ~exist(dirInfo,'dir')
    disp(['Creating ' dirInfo])
    mkdir(dirInfo);
end;


mi0=meCreateMicrographInfoStruct14;  % Get the mi file template
% Put in the invariant parameters
mi0.kV=lp.kV;
mi0.camera=lp.cameraIndex;
mi0.cpe=lp.cpe;
mi0.originalBasePath=AddSlash(rootPath);
mi0.basePath=mi0.originalBasePath;
mi0.moviePath='';
mi0.gainRefName='';
mi0.infoPath=dirInfo;
mi0.pixA=lp.defaultPixA;
    if mi0.kV>200
        mi0.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else
        mi0.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;

totalNFiles=0;
miNames={};

% Loop over subdirectories
% We pick up a list mNames of image names
for paIndex=1:numel(sdirName)
    dirImgs=sdirName{paIndex};
    pathImgs=AddSlash(dirImgs);
    disp(' ');
    disp(['Directory: ' micrDir pathImgs]);
    d=dir([micrDir pathImgs]);
    numFiles=0;
    sNames={};  % names in subdirectory
    for i=1:numel(d)  % scan the subdirectory and get image names
        [~,nm,ex]=fileparts(d(i).name);
        if any(strcmpi(ex,lp.imageExtension))
            numFiles=numFiles+1;
            sNames{numFiles}=d(i).name;
            sPaths{numFiles}=[micrDir pathImgs];
        end;
    end;
    
    % ----- Actual work is done here -----
    disp([num2str(numFiles) ' images to process']);
    for findex=1:lp.numExposures:numFiles
        mi=mi0;
        mi.imagePath=[micrDir pathImgs];
        % Construct the baseFilename
        [~,baseName,~]=fileparts(sNames{findex});
        %         Construct the mi filename here
        if lp.includeDirectoryName
            mi.baseFilename=[dirImgs '_' baseName];
        else
            mi.baseFilename=baseName;
        end;
        %         Read the image file header to get the pixel size, if possible
        if lp.checkImagePixA
            [~,s]=ReadEMFile([mi.imagePath sNames{findex}]);
            if s.pixA>0
                mi.pixA=s.pixA;
            elseif findex==1  % do this only once
                disp(['Setting mi.pixA to the default ' num2str(lp.defaultPixA)]);
            end;
        elseif findex==1  % do this only once
            disp(['Setting mi.pixA to the default ' num2str(lp.defaultPixA)]);
        end;
        
        %         Set the image filenames
        for ifocal=1:lp.numExposures
            mi.imageFilenames{ifocal}=sNames{findex+ifocal-1};
        end;
        %         Initialize the log.
        mi.log=cell(0,1);
        mi.log{end+1,1}=['CCDCreateInfoFiles ' TimeStamp];

%         Write the mi file
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

%         Accumulate the name list
        totalNFiles=totalNFiles+1;
        miNames{totalNFiles,1}=savedName;
        disp([str ': ' savedName]);
    end;
end;
disp([num2str(totalNFiles) ' files total.']);
disp(' ');

