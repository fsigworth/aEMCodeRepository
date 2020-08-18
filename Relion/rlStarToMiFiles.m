function rlStarToMiFiles(starName,pars)
% function mi=rlStarToMiFiles(starName,pars)
% function mi=rlStarToMiFiles(starCells,pars)
% From a micrographs_ctf.star file, read the image filenames and ctf parameters
% and create a set of mi files. We create the
% Info/ directory in the basePath (default the current directory)
% to contain the mi files.
% If the first argument is missing or is '' a file selector is put up.
% If the first argument starCells is a cell array, it contains the outputs from
% ReadStarFile() so the file doesn't have to be read again.
% The optinal second argument pars is a struct; see below for defaults.

if nargin<1
    starName=''; % we may have to put up the file selector.
end;
if nargin<2
    pars=struct; % use all defaults.
end;

dpars=struct; % Set up the defaults.

dpars.basePath=pwd; % assume we're in the relion project directory
dpars.cameraIndex=5; % K2
dpars.cpe=0.8;  % counts per electron, 0.8 for K2 counting mode, but
%  0.2 for superres image that is binned by MotionCor2.
% ! For Falcon3: cameraIndex=6, I think cpe=64.
dpars.dose=60; % Approx total movie dose in e/A^2. We have to guess this
% because MotionCor2 scaling doesn't allow the total dose to be calculated.
dpars.estimateStatsFromImage=0; % 1: don't use the above, estimate from image spectrum
dpars.nFrames=40;
dpars.BFactor=60; % Used by my CTF functions. Not critical.
dpars.changeImagePath=0; % change from the path given in the star file
dpars.imagePath='Micrographs/';
dpars.defaultImageSize=[0 0]; % Value to insert if we can't read the
% mrc file to get header information. If zeros, Vesicle_finding_GUI will
% assign this upon reading the image.
% dpars.readMicrographForScale=false; % If true, uses micrograph statistics
% the actual image. Slower because each file is read.
dpars.skipMissingMrcs=true; % Skip over any absent micrographs
dpars.writeMiFile=1; % Write out each mi.txt file.
dpars.setProcImage=0; % Set proc path to image path (if not writing merged images)
dpars.writeMergedImage=0;
dpars.writeMergedSmall=1;
dpars.writeJpeg=1;
dpars.dsSmall=4; % downscaling for Small and jpeg
pars=SetDefaultValues(dpars,pars,1); % 1 means check for undefined fieldnames.

cd(pars.basePath);

jpegPath='jpeg/';

if isa(starName,'cell') % Contains output from ReadStarFile()
    names=starName{1};
    dat=starName{2};
else % starName might be a string. If empty, we get starName from a file
    %     selector. We do not change our working directory.
    if numel(starName)<1
        disp('Getting a star file');
        [sName,sPath]=uigetfile('*.star');
        if isnumeric(sPath) % user clicked Cancel
            return
        end;
        starName=[sPath sName];
    end;
    
    [names,dat]=ReadStarFile(starName);
end;

[~,~,~,nLines]=rlStarLineToMi(names,dat,0);
disp([num2str(nLines) ' entries.']);
if pars.skipMissingMrcs
    disp('Skipping lines with no micrographs');
end;

%%
first=true;
skipCount=0;
for i=1:nLines
    [readOk,micFound,mi]=rlStarLineToMi(names,dat,i,pars);
    if ~readOk
        error(['Error reading star file data at line ' num2str(i)]);
    end;
    if first
        if pars.writeMiFile 
            CheckAndMakeDir(mi.infoPath,1);
            disp('Written:');
        end;
        if pars.writeMergedImage
            CheckAndMakeDir(mi.procPath,1);
        end
        if pars.writeMergedSmall % We now use Merged_sm as a new directory.
            CheckAndMakeDir(mi.procPath_sm,1);
        end;
        if pars.writeJpeg
            CheckAndMakeDir(jpegPath,1);
        end;
        first=false;
    end;
    if ~micFound && pars.skipMissingMrcs % silently skip these lines.
        skipCount=skipCount+1;
        continue;
    elseif pars.skipMissingMrcs
        if skipCount>0
            disp([num2str(skipCount) ' files skipped, up to ' mi.baseFilename]);
            skipCount=0;
        else % Estimate scaling
            mi=rlSetImageScale(mi,3,pars.nFrames);
        end;
    end;
    scaleMode=1+(pars.estimateStatsFromImage>0);
    if pars.writeMergedImage || pars.writeMergedSmall || pars.writeJpeg
        if pars.writeMergedImage
                    [mi,m,origSize]=rlSetImageScale(mi,scaleMode,pars.nFrames);
            micName=[mi.procPath mi.baseFilename 'm.mrc'];
            WriteMRC(m,mi.pixA,micName);
            padSize=mi.imageSize;
        else % we'll rely on reading the original micrograph
                    [mi,m,origSize]=rlSetImageScale(mi,scaleMode,pars.nFrames);
            padSize=mi.imageSize;
            mi.imageSize=origSize;
        end;
        if pars.writeMergedSmall || pars.dwriteJpeg
            smallName=[mi.procPath_sm mi.baseFilename 'ms.mrc'];
            jpegName=[jpegPath mi.baseFilename 'ms.jpg'];
                
            ms=Downsample(Crop(m,padSize),padSize/pars.dsSmall);
            if pars.writeMergedSmall
                WriteMRC(ms,mi.pixA*pars.dsSmall,smallName);
            end;
            if pars.writeJpeg
                WriteJpeg(ms,jpegName);
            end;
        else
            ms=BinImage(m,pars.dsSmall);
        end;
            imags(ms);
            title([num2str(i),': ' mi.baseFilename],'interpreter','none');
            drawnow;
    elseif  micFound
        [mi,~,origSize]=rlSetImageScale(mi,scaleMode,pars.nFrames);
        mi.imageSize=origSize;
    else         % Do dead reckoning (mode 3) if micrograph is not available.
        [mi,m,origSize]=rlSetImageScale(mi,3,pars.nFrames);
        mi.imageSize=origSize;
    end;
        miName=WriteMiFile(mi);
        disp([num2str(i) ': ' miName]);
end;

