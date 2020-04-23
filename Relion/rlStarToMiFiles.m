function rlStarToMiFiles(starName,pars)
% function mi=rlStarToMiFiles(starName,pars)
% function mi=rlStarToMiFiles(starCells,pars)
% From a micrographs_ctf.star file, read the image filenames and ctf parameters
% and create a set of mi files. We create the
% Info/ directory to contain the mi files.
% If the first argument is missing or is '' a file selector is put up.
% if pars. writeImage is true, Merged/ is created and contains padded and scaled
% micrograph files.
% If the first argument starCells is a cell array, it contains the outputs from
% ReadStarFile() so the file doesn't have to be read again.
if nargin<1
    starName=''; % we may have to put up the file selector.
end;
if nargin<2
    pars=struct; % use all defaults.
end;
dpars=struct;

dpars.basePath=pwd; % assume we're in the relion project directory
dpars.cameraIndex=5; % K2
dpars.cpe=0.2;  % counts per electron, 0.8 for K2 counting mode, but
%  0.2 for superres image that is binned by MotionCor2.
% ! For Falcon3: cameraIndex=6, I think cpe=64.


dpars.dose=60; % Approx total movie dose in e/A^2. We have to guess this
% because of an error in MotionCor2 scaling.
dpars.BFactor=60; % Used by my CTF functions. Not critical.
dpars.changeImagePath=1;
dpars.imagePath='';
dpars.readMicrographForScale=0; % Gets micrograph statistics directly from
% images, but greatly slows down execution.
dpars.skipIfNoImage=0; % if yes, we don't care if the image file is missing,
% we create the mi structure anyway.
dpars.writeMiFile=1; % Write it out.

dpars.writeFullSize=0; % write out full-size *.m image
dpars.writeDownsampled=0;
dpars.ds=8;  % downsampling factor for 'small' image

pars=SetDefaultValues(dpars,pars,1); % 1 means check for undefined fieldnames.

cd(pars.basePath);

if isa(starName,'cell')
    names=starName{1};
    dat=starName{2};
else % Get starName from a file selector. We do not change our path.
    if numel(starName)<1
        disp('Getting a star file');
        [starName,starPath]=uigetfile('*.star');
        if isnumeric(starPath) % user clicked Cancel
            return
        end;
        pars.starName=[starPath starName];
    end;
    [names,dat]=ReadStarFile(pars.starName);
end;

[~,~,~,nLines]=rlStarLineToMi(names,dat,0);
disp([num2str(nLines) ' entries.']);
if pars.readMicrographForScale
    disp('Reading micrographs to set scaling.');
end;
first=true;
for i=1:nLines
    [ok,mi,m]=rlStarLineToMi(names,dat,i,pars);
    if ok && pars.writeMiFile
        if first
            CheckAndMakeDir(mi.infoPath);
            disp('Written:');
            first=false;
        end;
        miName=WriteMiFile(mi);
        disp([num2str(i) ' ' miName]);
        if numel(m)>0 % we read the image, let's display it.
            nsz=NextNiceNumber(mi.imageSize);
            imags(BinImage(Crop(m,nsz,0,mean(m(:))),4));
            title(mi.baseFilename,'interpreter','none');
            drawnow;
        end;
    end;
end;

