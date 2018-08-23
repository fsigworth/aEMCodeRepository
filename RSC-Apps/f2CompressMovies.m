function f2CompressMovies(miNames,outRoot,compressMovies,compressImages)
% Read each mi file and compress the movie file(s), placing them into a
% parallel experiment tree starting with outRoot.

if nargin<3
    compressMovies=true;
end;
if nargin<4
    compressImages=false;
end;

moviePars=struct;
moviePars.snrRatio=100;  % target noise-to-noise ratio
moviePars.lfCutoff=.05; % radius of frequency domain not fitted

moviePars.displayOn=0;

imagePars=moviePars;
imagePars.lfCutoff=.1;

% For merged images:
% pars.snrRatio=100;  % target noise-to-noise ratio
% pars.lfCutoff=.2; % radius of frequency domain not fitted

if nargin<1 || numel(miNames)<1
    [miNames,pa]=uigetfile('*mi.*','Select mi files','multiselect','on');
    if isnumeric(miNames)
        return
    elseif ~isa(miNames,'cell')
        miNames={miNames};
    end;
    [rootPath,~]=ParsePath(pa);
    cd(rootPath);
end;
if nargin<2 || numel(outRoot)<1
    outRoot=uigetdir('','Select the output experiment directory');
    if isnumeric(outRoot)
        return
    end;
    moviePars.displayOn=(isdeployed || moviePars.displayOn) && ~useParFor;  % we put up the figure window if running alone.
end;
outRoot=AddSlash(outRoot);
nmi=numel(miNames);

% Create output directories and load mi files
CheckAndMakeDir([outRoot 'Info/'],1);
mis=cell(nmi,1);
for i=1:nmi
    disp([num2str(i) ': ' miNames{i}]);
    mi=ReadMiFile(miNames{i});
    mis{i}=mi;
    if compressMovies
        CheckAndMakeDir([outRoot mi.moviePath],1);
    end
    if compressImages
        CheckAndMakeDir([outRoot mi.imagePath],1);
    end;
end;

% Compression loop
tic
for i=1:nmi
    mi=mis{i};
    if exist([outRoot miNames{i}],'file') % mi file already exists in output path
        mi2=ReadMiFile([outRoot miNames{i}]);
    else
        mi2=mi;
    end;
    mi2.basePath=outRoot;
    if compressMovies
        movieNames=mi.movieFilename;
        if ~isa(movieNames,'cell')
            movieNames={movieNames};
        end;
        mi2.movieFilename=movieNames;  % allocate the cell array
        for j=1:numel(movieNames)
            finName=[mi.moviePath movieNames{j}];
            [~,nm,~]=fileparts(movieNames{j});
            foutName=[nm '.z.tif'];
            mi2.movieNames{j}=foutName;
            foutNameFull=[mi2.basePath mi2.moviePath foutName];
            disp(['Reading ' finName]);
            [m,s,fok]=ReadMovie(finName);  %% Note that we assume the movie file contains pixA.
            if fok
                disp(['Writing ' foutNameFull]);
                WriteZTiff(m,s.pixA,foutNameFull,moviePars);
            end;
        end;
    end;
    if compressImages
        imageNames=mi.imageFilenames;
        mi2.imageFilenames=imageNames;  % allocate the cell array
        for j=1:numel(imageNames)
            finName=[mi.imagePath imageNames{j}];
            [~,nm,~]=fileparts(imageNames{j});
            foutName=[nm '.z.tif'];
            mi2.imageFilenames{j}=foutName;
            foutNameFull=[mi2.basePath mi2.imagePath foutName];
            disp(['Reading ' finName]);
            [m,pixA,fok]=ReadEMFile(finName);  %% Note that we assume the image file contains pixA.
            if fok
                disp(['Writing ' foutNameFull]);
                WriteZTiff(m,pixA,foutNameFull,imagePars);
            end;
        end;
    end;
    outMiName=[mi2.basePath miNames{i}];
    WriteMiFile(mi2,outMiName);
end;
toc

