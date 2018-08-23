function miNames=spCreateInfoFiles(ciPars)
% function miNames=f2CreateInfoFiles(ciPars)
% The program creates a directory Info/ at the same level as the
% micrographs directory.  Into this we write an mi file corresponding to
% each micrograph.  The name is constructed from the characters baseNameChars
% of the micrograph, and the suffix 'mi'.
% Into each mi file are stored the relevant filenames and paths, and also
% values are set for (defaults are in parentheses)
% gp.camera     (5)  % K2
% gp.cpe        (0.8)
% gp.kV         (300)
% The names of the mi files are stored in the cell array gp.miNames.
% Typical values are
% gp.miNames{1}='Info/sq05_001_Jun10_00.01.8mi.txt'
% gp.miNames{2}='Info/sq05_002_Jun10_00.01.38mi.txt'

disp('spCreateInfoFiles');

% Set defaults for local parameters
lp.batchMode=0;  % batch mode is forced anyway if no gui
lp.overwrite=0;
lp.cameraIndex=5;  % K2 camera.
lp.cpe=0.8;        % counts/electron
lp.kV=300;
lp.defaultPixA=1.2156;
lp.doses=20;

imageExtension={'.mrc' '.tif'};
dirInfo='Info/';
baseNameChars=1:10;

% Overwrite these with values from the ciPars struct if present.
if nargin<1
    ciPars=struct;  % default is an empty struct
end;
lp=SetOptionValues(lp,ciPars);

% Have the user select the movie folder, or else find it.
disp('Getting the directory containing micrographs');
framesPath=uigetdir('.','Select a directory containing micrograph files');
if isnumeric(framesPath)
    return
end;
[rootPath,framesDir]=ParsePath(framesPath);
cd(rootPath);
%
% We pick up a list sNames of file names
disp(['Directory: ' framesDir]);
d=dir(framesDir);
numFiles=0;
sNames={};  % names in subdirectory
for i=1:numel(d)  % scan the directory and get movie names
    [~,nm,ex]=fileparts(d(i).name);
    if any(strcmpi(ex,imageExtension)) && nm(1)~='.' % check extension
        numFiles=numFiles+1;
        sNames{numFiles}=d(i).name;
    end;
end;
if numel(sNames)==0
    error(['No micrograph files found in ' framesDir]);
end;

if ~exist(dirInfo,'dir')
    disp(['Creating ' dirInfo])
    mkdir(dirInfo);
end;

%% ------actual work is done here-------

mi0=meCreateMicrographInfoStruct14;  % Get the mi file template
% Put in the invariant parameters
mi0.kV=lp.kV;
mi0.camera=lp.cameraIndex;
mi0.cpe=lp.cpe;
mi0.doses=lp.doses;
mi0.originalBasePath=AddSlash(rootPath);
mi0.basePath=mi0.originalBasePath;
mi0.infoPath=dirInfo;
mi0.moviePath='';
% Use the G&G damage model
if mi0.kV>200
    mi0.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
else
    mi0.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
end;

totalNFiles=0;
miNames={};


disp([num2str(numFiles) ' files to process']);
for fIndex=1:numFiles
    % Construct the output name
    %         serNumber=sprintf(lp.serNumberFormat,fileNumber);
    mi=mi0;
    % Construct the baseFilename
    [~,bName,ex]=fileparts(sNames{fIndex});
    outBaseName=[bName(baseNameChars)];
    
    %         Construct the mi filename here
    %         if outBaseName(end)=='4'  %% special case:
    %             outBaseName(end)=[];  %% remove trailing character for old packed mrc files.
    %         end;
    mi.baseFilename=outBaseName;
    miBasename=[dirInfo mi.baseFilename 'mi'];
    miTextName=[miBasename '.txt'];
    savedName=miTextName;  % default
    miExists=exist(miTextName,'file');
    if ~miExists || lp.overwrite % construct the new mi file
        
        mi.imagePath=AddSlash(framesDir);
        mi.imageFilenames(1)=sNames(fIndex);
        [m,pixA,ok]=ReadEMFile([framesDir sNames{fIndex}]);
        % Get the file information
        mi.imageSize=size(m);
        if ~isnan(pixA) && pixA>0.1
            mi.pixA=pixA;
        else
            mi.pixA=lp.defaultPixA;
        end;
        mi.doses=mean(m(:))/(mi.cpe*mi.pixA^2);
        mi.log={['spCreateInfoFiles ' TimeStamp]};
        
        %         Write the mi file.
        savedName=WriteMiFile(mi,1);
        totalNFiles=totalNFiles+1;
        if miExists
            str='overwritten';
        else
            str='written';
        end;
    else
        str='File exists';
    end;
    disp([str ': ' savedName]);
    miNames{fIndex,1}=savedName;
end; % for findex
disp([num2str(numFiles) ' micrographs total.']);
disp(' ');

