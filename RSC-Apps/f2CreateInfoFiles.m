function miNames=f2CreateInfoFiles(ciPars)
% function miNames=f2CreateInfoFiles(ciPars)
% Find movie files and create the corresponding mi (micrograph info) files,
% with parameters read from the optional ciPars structure.  The output is a
% cell array of the mi file names which have been written, including infoPath, e.g.
%       Info/sq0_0001_Jan14_23.59.59mi.txt
%
% Directory structure
%   The existing pair of movies (or sets, depending on ciPars.numExposures)
%   is assumed to reside in a directory like this:
% experiment/movie_frames/sq01_1/Sep24_18.52.48.mrc.
% In batch mode (either matlab was started with -nojvm, or
% ciPars.batchMode=1), we assume that we're in the experiment/ directory and
% this program looks for a directory whose name starts with 'movie'.  In
% interactive mode, a file selector is put up and you would select the
% 'movie_frames' directory.
% Either way, the subdirectories inside 'movie_frames' are searched for
% files with the extensions .mrc or .tif.  These are assumed to be movies.
% Thus in the experiment/ directory, we expect to find something like
%   movie_frames/sq05/Jun10_00.01.8.mrc  % first pair
%   movie_frames/sq05/Jun10_00.01.28.mrc
%               ...
%   movie_frames/sq08/Jun10_01.10.12.mrc % last pair
%   movie_frames/sq08/Jun10_01.10.31.mrc
%               ...
% The program creates a directory Info/ at the same level as the
% movie_frames directory.  Into this we write an mi file corresponding to
% each movie.  The name is constructed from the subdirectory name, a serial
% number, the first original movie name, and the suffix 'mi'.  In the example we will have
%   Info/sq05_001a_Jun10_00.01.8mi.txt
%   Info/sq05_001b_Jun10_00.01.38mi.txt
%               ...
%   Info/sq08_010a_Jun10_01.10.1mi.txt
%   Info/sq08_010b_Jun10_01.10.31mi.txt
%               ...
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
lp.numExposures=2;
lp.batchMode=0;  % batch mode is forced anyway if no gui
lp.overwrite=0;
lp.cameraIndex=6;  % 6 means Falcon2 camera.
lp.cpe=50;        % counts/electron
lp.kV=200;
lp.defaultPixA=1.247;
lp.checkMoviePixA=0;  % Don't read movie files to look for pixA.
lp.serNumberFormat='%04d';
movieExtension={'.mrc' '.tif'};
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
    framesPath=uigetdir('.','Select a directory containing directories of movie files');
    if isnumeric(framesPath)
        return
    end;
    [rootPath,framesDir]=ParsePath(framesPath);
    cd(rootPath);
else
    disp('Running batch mode.');
    d=dir;
    framesDir=[];
    for i=1:numel(d)
        if strncmpi(d(i).name,'movie',5) && d(i).isdir
            framesDir=AddSlash(d(i).name);
            break
        end;
    end;
    if numel(framesDir)<1
        error(['No movie directory in current dir: ' pwd]);
    end;
    rootPath=AddSlash(pwd);
end;

d=dir(framesDir);  % scan for subdirectories of movie_frames/
nd=0;
pathName={};
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
mi0.camera=lp.cameraIndex;
mi0.cpe=lp.cpe;
mi0.originalBasePath=AddSlash(rootPath);
mi0.basePath=mi0.originalBasePath;
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
for paIndex=1:numel(pathName)
    dirMovie=pathName{paIndex};
    disp(' ');
    disp(['Directory: ' framesDir dirMovie]);
    d=dir([framesDir dirMovie]);
    numFiles=0;
    sNames={};  % names in subdirectory
    for i=1:numel(d)  % scan the subdirectory and get movie names
        [~,nm,ex]=fileparts(d(i).name);
% disp(d(i).name);
% disp(ex);
        if any(strcmpi(ex,movieExtension)) && nm(1)~='.' % check extension
            numFiles=numFiles+1;
            sNames{numFiles}=d(i).name;
        end;
    end;
    
    disp([num2str(numFiles) ' files to process']);
    fileNumber=1;
    for fIndex=1:lp.numExposures:numFiles
        % Construct the output name
        serNumber=sprintf(lp.serNumberFormat,fileNumber);
        mi=mi0;
        % Construct the baseFilename
        [~,baseName,ex]=fileparts(sNames{fIndex});
        
        if (strcmp(ex,'.tif') && baseName(end)=='z')
            baseName(end)=[];  % remove the terminal 'z'
        end;
        
        if dirMovie(end)=='/'
            dirMovie(end)=[];  % prepare to concatenate the path name
        end;
        %         Construct the mi filename here
        outBaseName=[dirMovie '_' serNumber '_' baseName];
%         if outBaseName(end)=='4'  %% special case:
%             outBaseName(end)=[];  %% remove trailing character for old packed mrc files.
%         end;
        mi.baseFilename=outBaseName;
            miBasename=[dirInfo mi.baseFilename 'mi'];
            miTextName=[miBasename '.txt'];
            savedName=miTextName;  % default
            miExists=exist(miTextName,'file');
        if ~miExists || lp.overwrite % construct the new mi file
            
            mi.moviePath=[AddSlash(framesDir) AddSlash(dirMovie)];
            % Get the file information
            [~,s,ok]=ReadMovie([mi.moviePath sNames{fIndex}],1,0);
            if ok  % We got a valid file
                fileNumber=fileNumber+1;  % increment it only if a valid file.
                mi.pixA=s.pixA;
                mi.frameSets=[1 s.nz];
                mi.movieFilename{1}=sNames{fIndex};
                lastEIndex=min(lp.numExposures,numel(sNames)-fIndex+1);
                %                 Get up to numExposures of additional movies
                for eIndex=2:lastEIndex
                    fname=sNames{fIndex+eIndex-1};
                    [~,~,mok]=ReadMovie([mi.moviePath fname],1,0);
                    if ~mok
                        break;
                    end;
                    mi.movieFilename{eIndex}=fname;
                    mi.frameSets(eIndex,:)=[1 s.nz];
                end;
                mi.imageSize=[s.nx s.ny];
                %         Read the movie file header to get the pixel size, if possible

                %         Initialize the log.
                mi.log={['f2CreateInfoFiles ' TimeStamp]};

                %         Write the mi file.
                savedName=WriteMiFile(mi,1);
                totalNFiles=totalNFiles+1;
                miNames{totalNFiles,1}=savedName;
                if miExists
                    str='overwritten';
                else
                    str='written';
                end;
            else
                disp(['Not a valid micrograph: ' mi.moviePath sNames{fIndex}]);
                str='Not written';
            end; % if ok
        else
            str='Already exists';
        end;  % if exists
        disp([str ': ' savedName]);
    end; % for findex
end; % for paIndex
disp([num2str(totalNFiles) ' movie sets total.']);
disp(' ');

