%function miNames=k2CreateInfoFiles2(ciPars)
ciPars=struct;
%ciPars.defaultPixA=1.533;
ciPars.defaultPixA=1.333;
ciPars.readOnly=0;
ciPars.overwrite=0;
ciPars.getAllNames=0;
ciPars.numExposures=2;
ciPars.kV=300;
ciPars.batchMode=0;
ciPars.displayOn=0;
% ciPars.rawDir='/ysm-gpfs/pi/cryoem/krios/170808/KvLipo/';
ciPars.rawDir='~/EMWork/Hideki/170808/KvLipo/';

nargin=1;

% modified to allow rawDir distinct from the experiment path fs 20Aug17


%function miNames=k2CreateInfoFiles(ciPars)
% Find movie files and create the corresponding mi (micrograph info) files,
% with parameters read from the optional ciPars structure.  The output is a
% cell array of the mi file names.
%
% Directory structure
%   The existing movie is assumed to reside in a directory like this:
% experiment/movie_frames/sq01_1/Sep24_18.52.48.mrc (or else .tif
% extension').  In batch mode (either matlab was started with -nojvm, or
% gp.batchMode=1), we assume that we're in the experiment/ directory and
% this program looks for a directory whose name starts with 'movie'.  In
% interactive mode, a file selector is put up and you would select the
% 'movie_frames' directory.
% Either way, the subdirectories inside 'movie_frames' are searched for
% files with the extensions .mrc or .tif.  These are assumed to be movies.
% Thus in the experiment/ directory, we expect to find something like
%   movie_frames/sq05/Jun10_00.01.8.tif
%   movie_frames/sq05/Jun10_00.01.38.tif
%               ...
%   movie_frames/sq08/Jun10_01.10.1.tif
%   movie_frames/sq08/Jun10_01.10.31.tif
%               ...
% In interactive mode we ask for a gain reference file with extension .dm4
% In batch mode we search all of the directories sq01/, sq02/ etc. until we
% find a .dm4 file.  This is taken to be the gain reference.
%
% The program creates a directory Info/ at the same level as the
% movie_frames directory.  Into this we write an mi file corresponding to
% each movie.  The name is constructed from the subdirectory name, a serial
% number, the original movie name, and the suffix 'mi'.  In the example we will have
%   Info/sq05_001_Jun10_00.01.8mi.txt
%   Info/sq05_002_Jun10_00.01.38mi.txt
%               ...
%   Info/sq08_001_Jun10_01.10.1mi.txt
%   Info/sq08_002_Jun10_01.10.31mi.txt
%               ...
% Into each mi file are stored the relevant filenames and paths, and also
% values are set for (defaults are in parentheses)
% gp.camera     (5)
% gp.cpe        (0.8)
% gp.kV         (200)
% gp.defaultPixA       (1.247)
%   Note that the defaultPixA is used only in the case (generic .tif) that
%   no pixA value is stored.
% The names of the mi files are stored in the cell array gp.miNames.
% Typical values are
% gp.miNames{1}='Info/sq05_001_Jun10_00.01.8mi.txt'
% gp.miNames{2}='Info/sq05_002_Jun10_00.01.38mi.txt'


% % disp('k2CreateInfoFiles');

% Set defaults for local parameters
lp.batchMode=0;
lp.overwrite=0;
lp.getAllNames=0;  % return not just created files, but existing file names too.
lp.readOnly=0;     % don't create new files.
lp.cameraIndex=5;  % 5 means K2 camera.
lp.cpe=0.8;        % counts/electron, here reflects the counting efficiency
lp.kV=200;
lp.defaultPixA=1.533;
lp.checkMoviePixA=0;  % Don't read movie files to look for pixA.
lp.numExposures=2;  % single movie with a defocus jump in the middle.
lp.indexFormat='%04d';
lp.displayOn=1;
lp.rawDir=''; % where to find the framesDir.  Default is in the experiment directory.

%movieExtension={'.tif' '.mrc'};
movieExtension={'.tif'};
refExtension='.dm4';
dirInfo='Info/';

% Overwrite these with values from the ciPars struct if present.
if nargin<1
    ciPars=struct;  % default is an empty struct
end;
lp=SetOptionValues(lp,ciPars);

lp.batchMode=lp.batchMode || ~usejava('desktop');  % Matlab's gui is not active, must use batch mode
% lp.batchMode=0;  %%%%%

% Have the user select the directories, or else find the raw data directory.
if ~lp.batchMode  % Matlab's GUI is active
     disp('Getting the experiment directory');
    lp.rootPath=uigetdir('../','Select the experiment directory');
    if isnumeric(lp.rootPath)
        return
    end;
    cd(lp.rootPath);
    disp('Getting the movie_frames directory');
    framesDir=AddSlash(uigetdir('.','Select the movie_frames directory'));
    if isnumeric(framesDir)
        return
    end;
else
    % framesDir is either lp.rawDir, or [lp.rawDir 'movie_frames']
    % movies are either in framesDir or in subdirectories under framesDir.
    % Hence we handle any of these file structures:
    % 170808/*.tif
    % 170808/movie_frames/*.tif
    % 170808/movie_frames/sq0*/*.tif (can be multiple sq directories)
    % here framesDir='170808/movie_frames/' and pathName{1}='sq01/'
    d=MyDir(lp.rawDir);
    framesDir=AddSlash(lp.rawDir); % fefault for first case of the above
    for i=1:numel(d) % search for movie* as a directory name; use it if found.
        if strncmpi(d(i).name,'movie',5) && d(i).isdir
            framesDir=AddSlash(d(i).name);
            break
        end;
    end;
%     if numel(framesDir)<1
%         error('Couldn''t find a movie directory');
%     end;
%     rootPath=AddSlash(pwd);
end;
%%
% Scan for subdirectories
        d=MyDir(framesDir);
nd=0;
pathName={''};  % default is no subdirectories
for i=1:numel(d)
    if d(i).isdir && d(i).name(1)~='.'
        nd=nd+1;
        pathName{nd}=AddSlash(d(i).name);
    end;
end;

if nd>0  % subdirectories found
%     disp([num2str(nd) ' directories found']);
else
    disp(['No enclosed directories found at ' lp.rawDir framesDir]);
    disp('Looking for movie files there.');
end;

%% Pick up a gain reference which will be used for all movies

if ~lp.batchMode
    oldDir=pwd;
    cd(framesDir);
    disp('Selecting the gain reference');
    [refName, refPath]=uigetfile('*.dm4','Select a gain reference');
    cd(oldDir);
    if ischar(refName)
        %         refPath is returned as an absolute path.  Make it relative
        refPath=AddSlash(GetRelativePath(refPath,pwd));
        refName=[refPath refName];
    end;
else  % search for the first .dm4 file in the movies folder
    % Loop over subdirectories and find the first reference.
    refName=0;  % default is no reference file.
    paIndex=1;
    while isa(refName,'numeric') && paIndex<=numel(pathName)
        dirMovie=pathName{paIndex};
        refPath=[framesDir dirMovie];
        disp(['Directory: ' refPath]);
        d=MyDir(refPath);
        for i=1:numel(d)
            [~,nm,ex]=fileparts(d(i).name);
            if any(strcmpi(ex,refExtension)) && nm(1)~='.' % we take .dm4 files
                refName=[refPath d(i).name];
                break;
            end;
        end;
        paIndex=paIndex+1;
    end;
end;

if ischar(refName) % user selected something
    disp(['Assigning the initial gain reference: ' refName]);
else
    disp('No gain reference used.');
    refName='';
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
mi0.originalBasePath=AddSlash(pwd);
mi0.basePath=mi0.originalBasePath;
mi0.gainRefName=refName;
mi0.infoPath=dirInfo;
mi0.pixA=lp.defaultPixA;
% Use the G&G damage model
    if mi0.kV>200
        mi0.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else
        mi0.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;
  
totalNFiles=0;
miNames={};

% Loop over subdirectories
% We pick up a list mNames of movie names and also the latest reference filename
dirs=cell(numel(pathName),1);
for paIndex=1:numel(pathName)
    dirMovie=pathName{paIndex};
    disp(' ');
    disp(['Directory: ' framesDir dirMovie]);
    d=MyDir([framesDir dirMovie]);
    dirs{paIndex}=d;
    numFiles=0;
    sNames={};  % names in subdirectory
    for i=1:numel(d)  % scan the subdirectory and get movie names
        [~,nm,ex]=fileparts(d(i).name);
        if any(strcmpi(ex,movieExtension)) && nm(1)~='.' % we take .mrc, .tif files
            numFiles=numFiles+1;
            sNames{numFiles}=d(i).name;
        end;
    end;
    
    numMis=floor(numFiles/lp.numExposures);
    disp([num2str(numFiles) ' movie files total.']);
    disp([num2str(numMis) ' mi files total.']);
    
    for miIndex=1:numMis;
        fIndex=(miIndex-1)*lp.numExposures+1; % pick up first of exposure pair

        mi=mi0;
        % Construct the baseFilename
        [~,baseName,~]=fileparts(sNames{fIndex});
        % Construct the output name
        serNumber=sprintf(lp.indexFormat,miIndex);
        if dirMovie(end)=='/'  % un-do the AddSlash operation.
            dirMovie(end)=[];
        end;
        %         Construct the mi filename here
        outBaseName=[dirMovie '_' serNumber '_' baseName];
        %         if outBaseName(end)=='4'  %% special case:
        %             outBaseName(end)=[];  %% remove trailing character for old packed mrc files.
        %         end;
        mi.baseFilename=outBaseName;
        mi.moviePath=[AddSlash(framesDir) AddSlash(dirMovie)];
        %         Read the movie file header to get the pixel size, if possible
        mi.frameSets(1,:)=[1 inf];
        if lp.checkMoviePixA
            [~,s,mok]=ReadMovie([mi.moviePath sNames{fIndex}],1,0);
            if isfield(s,'pixA') && s.pixA>0
                mi.pixA=s.pixA;
            elseif fIndex==1  % do this only once
                disp(['Setting mi.pixA to the default ' num2str(lp.defaultPixA)]);
            end;
        elseif fIndex==1  % do this only once
            disp(['Setting mi.pixA to the default ' num2str(lp.defaultPixA)]);
        end;
        
        
        miBasename=[dirInfo mi.baseFilename 'mi'];
        miTextName=[miBasename '.txt'];
        savedName=miTextName;  % default
        str='';
        miExists=exist(miTextName,'file');  % don't already have a mi.txt file
%         if ~miExists  % check for an mi.mat file
%             miMatName=[miBasename, '.mat'];
%             if exist(miMatName,'file')
%                 mi=ReadMiFile(miMatName);
%                WriteMiFile(mi,miTextName);  % convert an existing .mat to .txt
%                 str='Converted to text';
%                 miExists=true;
%             end;
%         else
%             str='Already exists';
%         end;
% 
        if miExists
            str='Exists. ';
        end;
        
        doWrite=(~miExists || lp.overwrite) && ~lp.readOnly;
        doEnter=doWrite || lp.getAllNames;

        
        %         Enter the movie filename(s)
        if lp.numExposures==1
            mi.movieFilename=sNames{fIndex};   % a simple string
        else
            mi.movieFilename{1}=sNames{fIndex};  % first of a cell array
            for i=2:lp.numExposures
                if numel(sNames)>=fIndex+i-1
                    fname=sNames{fIndex+i-1};
                    if doWrite
                        [~,~,mok]=ReadMovie([mi.moviePath fname],1,0);
                    else
                        mok=true;
                    end;
                    
                    if mok
                        mi.movieFilename(i)=sNames(fIndex+i-1);
                        mi.frameSets(i,:)=[1 inf];
                    else
                        doWrite=false;
                        doEnter=false;
                    end;
                else
                    doWrite=false;
                    doEnter=false;
                end;
            end;
        end;
        
        %         Initialize the log.
        mi.log=cell(0,1);
        mi.log{end+1,1}=['k2CreateInfoFiles ' TimeStamp];
        
        %         Write the mi file, if desired.

        % Write or make entries into the miNames array.
        if doWrite
            savedName=WriteMiFile(mi,1);
            if lp.overwrite && miExists
                str=[str '; overwritten'];
 %               str='overwritten.';
            else
                 str='Written.';
            end;
        else
            str='Not written.';
        end;
        
        if doEnter
            totalNFiles=totalNFiles+1;
            miNames{totalNFiles,1}=savedName;
            if lp.displayOn
            disp([str ' Entered: ' savedName]);
            end;
        else
            if lp.displayOn
            disp([str ' Not entered:' savedName]);
            end;
        end
     end;
end;
% if totalNFiles>0
    disp([num2str(totalNFiles) ' files total.']);
    disp(' ');
% end;


function list=MyDir(path)
if nargin<1 || numel(path)<1
    list=dir;
else
    list=dir(path);
end;
end
