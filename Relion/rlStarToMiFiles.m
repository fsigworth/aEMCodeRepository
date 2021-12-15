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
% dpars.cpe=0.8;  % counts per electron, 0.8 for K2 counting mode, but
dpars.cpe=14/4; % K3 is scaled up by 16, but mc2 handles superres incorrecly

%  0.2 for superres image that is binned by MotionCor2.
% ! For Falcon3: cameraIndex=6, I think cpe=64.
dpars.dose=50; % Approx total movie dose in e/A^2. We have to guess this
% because MotionCor2 scaling doesn't allow the total dose to be calculated.
dpars.estimateStatsFromNoise=0; % 1: don't use the above, estimate from image spectrum%
dpars.nFrames=34;
dpars.motionCorFrames=dpars.nFrames; % either 1, or the number of frames if using MotionCor2
dpars.scaleMode=0; % 0: K3 micrographs, already normalized; imageNormScale set to 1.
% 1: read micrograph and scale for k2
% 2: estimate stats from noise
dpars.BFactor=40; % Used by my CTF functions. Not critical.
dpars.noDamageModel=1; % No damage filtering.

dpars.changeImagePath=0; % change from the path given in the star file
dpars.imagePath='Micrographs/'; % ...new image path
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
dpars.writeJpegInv=1;  % Make -1 to reverse the contrast.
dpars.compFraction=0.3;
dpars.dsSmall=4; % downscaling for Small and jpeg
dpars.disHP=0;
dpars.disFc=.4;
dpars.lastLine=inf;
dpars.firstPeakAmp=.5;

pars=SetDefaultValues(dpars,pars,1); % 1 means check for undefined fieldnames.

cd(pars.basePath);

jpegPath='Jpeg/';
jpegInvPath='JpegInv/';

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
    
    [names,dat]=ReadStarFile(starName,1,100+pars.lastLine);
end;

[~,~,~,nLines]=rlStarLineToMi(names,dat,0);
disp([num2str(nLines) ' entries.']);
nLines=min(nLines,pars.lastLine);
disp([num2str(nLines) ' to process.']);

if pars.skipMissingMrcs
    disp('Skipping lines with no micrographs');
end;

%%
first=true;
skipCount=0;
oldImageName='';
for i=1:nLines
    newName=dat{end}.rlnMicrographName{i};
    if strcmp(newName,oldImageName)
        continue;
    else
        oldImageName=newName;
    end;
    disp(newName);
    [readOk,micFound,mi]=rlStarLineToMi(names,dat,i,pars);
    if ~readOk
        error(['Error reading star file data at line ' num2str(i)]);
    end;
    if first
        if pars.writeMiFile && skipCount==0
            CheckAndMakeDir(mi.infoPath,1);
            disp('Written...');
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
        if pars.writeJpegInv
            CheckAndMakeDir(jpegInvPath,1);
        end;
    end;
    if pars.skipMissingMrcs % silently skip these lines.
        if ~micFound
            skipCount=skipCount+1;
            continue;
        elseif skipCount>0
            disp([num2str(skipCount) ' files skipped, up to ' mi.baseFilename]);
            skipCount=0;
        end;
    end;
    
    %     scaleMode=1+(pars.estimateStatsFromNoise>0);
    writeSomething=pars.writeMergedImage || pars.writeMergedSmall || ...
        pars.writeJpeg || pars.writeJpegInv;
    if writeSomething
        [mi,m]=rlSetImageScale(mi,pars.scaleMode,pars.motionCorFrames);
        if pars.writeMergedImage
            micName=[mi.procPath mi.baseFilename 'm.mrc'];
            WriteMRC(m,mi.pixA,micName);
        end;
        ms=Downsample(Crop(m,mi.padImageSize),mi.padImageSize/pars.dsSmall);
        msDis=imscale(GaussFilt(ms,pars.disFc),256,1e-4);
        
        if pars.writeMergedSmall || pars.writeJpeg || pars.writeJpegInv
            smallName=[mi.procPath_sm mi.baseFilename 'ms.mrc'];
            jpegName=[jpegPath mi.baseFilename 'ms.jpg'];
            jpegInvName=[jpegInvPath mi.baseFilename 'msinv.jpg'];
            
            if pars.writeMergedSmall
                WriteMRC(ms,mi.pixA*pars.dsSmall,smallName);
            end;
            if pars.writeJpeg
                WriteJpeg(msDis,jpegName);
            end;
            if pars.writeJpegInv ~= 0
                ipars=struct;
                ipars.fHP=pars.disHP;
                ipars.firstPeakAmp=pars.firstPeakAmp;
                msInv=pars.writeJpegInv*rspCTFInverseFilter(ms,mi,pars.compFraction,ipars);
                msInvDis=GaussFilt(msInv,pars.disFc);
                msDis=WriteJpeg(msInvDis,jpegInvName);
            end;                
        end;
        imaga(msDis);
        title([num2str(i),': ' mi.baseFilename],'interpreter','none');
        drawnow;
    else % not writing anything, don't need to read the original image
        mi=rlSetImageScale(mi,pars.scaleMode,pars.motionCorFrames);
    end; % writeSomething
    %     else         % Do dead reckoning (mode 3) if micrograph is not available.
    %         [mi,m,origSize]=rlSetImageScale(mi,3,pars.motionCorFrames);
    %         mi.imageSize=origSize;
    %     end;
    if first
        disp(['The first image median, normScale: ' num2str([mi.imageMedian mi.imageNormScale])]);
    end;
    
    miName=WriteMiFile(mi);
    disp([num2str(i) ': ' miName]);
    first=false;
end;

