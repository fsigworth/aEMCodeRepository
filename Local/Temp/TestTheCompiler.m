function TestTheCompiler(ciPars)

% Set defaults for local parameters
lp.batchMode=0;
lp.overwrite=0;
lp.cameraIndex=5;  % 5 means K2 camera.
lp.cpe=0.8;        % counts/electron, here reflects the counting efficiency
lp.kV=200;
lp.defaultPixA=1.247;
lp.checkMoviePixA=0;  % Don't read movie files to look for pixA.

movieExtension={'.tif' '.mrc'};
refExtension='.dm4';
dirInfo='Info/';

% Overwrite these with values from the ciPars struct if present.
if nargin<1
    ciPars=struct;  % default is an empty struct
end;
lp=SetOptionValues(lp,ciPars);

lp.batchMode=lp.batchMode || ~usejava('desktop');  % Matlab's gui is not active, must use batch mode
lp.batchMode=0;  %%%%%

% Have the user select the movie folder, or else find it.
if ~lp.batchMode  % Matlab's GUI is active
% %     disp('Getting the framesPath');
    framesPath=uigetdir('.','Select a directory containing directories of movie files');
    if isnumeric(framesPath)
        return
    end;
    [rootPath,framesDir]=ParsePath(framesPath);
    cd(rootPath);
else
    d=dir;
    framesDir=[];
    for i=1:numel(d)
        if strncmpi(d(i).name,'movie',5) && d(i).isdir
            framesDir=AddSlash(d(i).name);
            break
        end;
    end;
    if numel(framesDir)<1
        error('Couldn''t find a movie directory');
    end;
    rootPath=AddSlash(pwd);
end;

figure;
text(.5,.5,rootPath);
