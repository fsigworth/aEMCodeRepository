% rtOnlineDisplay.m

batchMode=1;
inWorkingDir=1; % We are starting in the experiment directory
% motionCorMode='ReadStar'; % ReadStar = read star file created by Relion's motioncorr
motionCorMode='RunMC2';
writeNormalizedImages=1; % Normalize and write into the Merged/ directory.
graphicsOn=0;
dss=4; % downsample factor for ms.mrc files
skipExistingMis=1; % if mi file exists, skip processing.
maxTotalDose=0; % if nonzero, skips processing if > maxTotalDose

tempDir='Logs/';
imageDir='Micrograph/';
jpegDir='Jpeg/';
infoDir='Info/';
procDir='Merged/';
ourRefName='CountRefLocal.mrc';
gainRefRot=0;
kV=200;
searchMode=2; % 1: operate on latest; 2 start from beginning;
%               3: same as 2 but repeating
pixA=0.86;
cpe=6;
defaultFrameDose=1.25;

cameraType=7; % 5 is k2; 7 is Falcon3 in counting mode.
times=0;
relTimes=0;
defoci=0;
maxRepeat=100; % 500 seconds max. waiting before quitting.
maxRepeat=1;

mc2Pars=struct;
mc2Pars.throw=1;   % for MC2: default is zero.
mc2Pars.trunc=0;  % default is zero
mc2Pars.gainRefRot=0;
mc2Pars.patches=[1 1];
mc2Pars.defaultFrameDose=defaultFrameDose;
mc2Pars.tempDir=tempDir;

if batchMode
    inPath='Micrographs2/';
    namePattern='.tif';
    ourRefName='';
else
    % Use file selectors to set
    %   inPath (path to movies), namePattern (e.g. '.tif'), ourRefName ('' if absent)
    disp('Getting the first movie');
    [firstName,inPath]=uigetfile('*.*','Select a movie file');
    
    if ~inWorkingDir
        cd(inPath);  % Might need to use a base path
        inPath='./';
    else
        inPath=GetRelativePath(inPath);
    end;
    if numel(inPath)<2
        inPath='./';
    end;
    [pa,nm,ex]=fileparts(firstName);
    namePattern=ex
    %namePattern='.mrc'
    
    %
    % Look for a gain reference of the form *ref*.dm4
    d=dir(inPath); % same as movie dir
    found=0;
    for i=3:numel(d)
        [~,rnm,rex]=fileparts(d(i).name);
        if strcmp(rex,'.dm4') && numel(strfind(lower(rnm),'ref'))
            refName=d(i).name;
            refPath=GetRelativePath(inPath);
            found=i;
            break;
        end;
    end;
    if found
        gainRefName=[refPath refName];
    else
        disp('Getting the reference');
        [refName,refPath]=uigetfile('*.dm4','Select the reference');
        if isnumeric(refName) % nothing selected
            gainRefName='';
        else
            gainRefName=[GetRelativePath(refPath) refName];
        end;
    end;
    if numel(gainRefName)>1
        ref=ReadEMFile(gainRefName);
        WriteMRC(ref,0,ourRefName);
    end;
    
    if numel(gainRefName)>1
        disp(['Gain reference: ' gainRefName ' copied to ' ourRefName]);
    else
        disp('No gain reference');
        ourRefName='';
    end;
    
end
%%

CheckAndMakeDir(tempDir, 1);
CheckAndMakeDir(imageDir, 1);
CheckAndMakeDir(jpegDir, 1);
CheckAndMakeDir(infoDir,1);
CheckAndMakeDir(procDir,1);

mi0=meCreateMicrographInfoStruct14;
mi0.pixA=pixA;
% mi0.movieFilename{1}=nameList{end};
mi0.tempPath=tempDir;
mi0.imagePath=imageDir;
mi0.procPath=procDir;
mi0.jpegPath=jpegDir;
mi0.gainRefName=ourRefName;
mi0.gainRefRot=gainRefRot;
mi0.kV=kV;
mi0.cpe=cpe;
mi0.camera=cameraType;
mi0.moviePath='';
    if mi0.kV>200
        mi0.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else
        mi0.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;

% ---------------------Main loop---------------
doRepeat=1;
miIndex=0;
% fileIndex=0;
mov=struct;
mov.inPath=inPath;
mov.pattern=namePattern;
    mov.freeList=logical([]);
    mov.nameList={};
    mov.timeList=[];
    mov.fileIndex=0;

while doRepeat
    ok=0;
    repeatCount=0;
    while ~ok && repeatCount<maxRepeat
        mov=rtFindNextMovie(mov,searchMode);
        if ~mov.foundAFile
            if searchMode==2 % only once through
                disp('done.');
                return
            end;
            sprintf('.');
            pause(5);
        end;
        repeatCount=repeatCount+1;
    end;
    doRepeat=repeatCount<=maxRepeat;
    
    miIndex=miIndex+1;
    times(miIndex)=mov.timeList(mov.fileIndex);

    mi=mi0;
    mi.movieFilename{1}=mov.nameList{mov.fileIndex};
    [pa,nm,ex]=fileparts(mi.movieFilename{1});
    mi.baseFilename=nm;
    
    switch motionCorMode
        case 'ReadStar' % Read a star file created by Relion, in the same dir as the images.
            mi.imagePath=inPath; % we're reading images
            mi.imageFilenames{1}=[nm '.mrc'];  % new filename, if we're running MotionCor

            starName=[inPath mi.baseFilename '.star'];
            disp(starName);
            if ~exist(starName,'file')
                error('...file not found.');
            end;
            [dn,s]=ReadStarFile(starName);
            % we assume that s{1} is data_general; s{2} is data_global_shift
            mi.imageSize=[s{1}.rlnImageSizeX s{1}.rlnImageSizeY];
            mi.pixA=s{1}.rlnMicrographOriginalPixelSize*s{1}.rlnMicrographBinning;
            mi.kV=s{1}.rlnVoltage;
            
            nl=numel(s{2}.rlnMicrographShiftX);
            mi.frameDose=ones(nl,1)*defaultFrameDose;
            mi.frameSets=[s{1}.rlnMicrographStartFrame nl];
            mi.doses=sum(mi.frameDose);
            mi.frameShifts=cell(1,1);
            mi.frameShifts{1}=single([s{2}.rlnMicrographShiftX s{2}.rlnMicrographShiftY]);
            
            microName=[mi.imagePath mi.baseFilename '.mrc'];
            [m,sh0]=ReadMRC(microName);
            %%%%%        WriteMRC(m,mi.pixA,[mi.imagePath mi.imageFilenames{1}]);
            %        return
        otherwise
            % Skip if info file already exists
            miName=[infoDir mi.baseFilename 'mi.txt'];
            if exist(miName,'file')
                if skipExistingMis
                    disp([miName ' exists. Skipped.']);
                    continue;
                end;
                mi1=ReadMiFile(miName);
                if maxTotalDose
                if all(mi1.doses< maxTotalDose) % reasonable dose
                    disp(['dose is ' num2str(mi1.doses) ', skipping']);
                        continue;
                else
                    disp(['dose is ' num2str(mi1.doses) ', reprocessing.']);
                end;
                end;
            end;
            mi.imageFilenames{1}=[nm 'ali.mrc'];  % new filename, if we're running MotionCor
            mi.moviePath=inPath;
            % ------Run MC2 here------
            mi=rtMC2Runner(mi,mc2Pars);
            if isnumeric(mi) % MC2 run was unsuccessful.
                disp('MotionCor2 run was unsuccessful.');
                %        return;
                break;
            end;
            m=ReadMRC([mi.imagePath mi.imageFilenames{1}]);
    end;
    n0=size(m);
    n=NextNiceNumber(n0,5,8);
    m1=Crop(m,n,0,mean(m(:))); % expand to nice size
    m1s=Downsample(m1,n/dss);
    if writeNormalizedImages
        % Relion's normalization has correct DC component
        me1=mean(m1(:));
        m1=m1/me1-1;
        m1s=m1s/me1-1;
        mergeNm=[mi.procPath mi.baseFilename];
        disp(['Writing ' mergeNm ' m.mrc and ms.mrc']);
        WriteMRC(m1,mi.pixA,[mergeNm 'm.mrc']);
        WriteMRC(m1s,mi.pixA*dss,[mergeNm 'ms.mrc']);
    end;
    
    %%%%%
    
    if graphicsOn
        mysubplot(122);
        imags(m1s);
        axis off;
        title(mi.movieFilename{1}, 'interpreter','none');
        
        mysubplot(2,4,1);
        plot(mi.frameShifts{1});
        grid on;
        ylabel('Shift, pixels');
        drawnow;
        
    else
        s=struct;
        s.exec='mysubplot(122)';
        s.image=m1s;
        s.axis='off';
        s.title=mi.movieFilename{1};
        disDat.mainImage=s;
        
        s=struct;
        s.exec='mysubplot(2,4,1)';
        s.ploty=mi.frameShifts{1};
        %     s.exec('grid on');
        %--    s.ylabel('Shift, pixels');
        disDat.frameShifts=s;
    end;
    
%
%-----------run Gctf----------------
    gPars=struct;
    gPars.tempDir=tempDir;
    [mi,epaVals,ctfImage,ctfVals]=rtGctfRunner(mi,gPars);
    if isnumeric(mi)
        disp('Skipping');
        break
    else
        disp(['defocus ' num2str(mi.ctf(1).defocus)]);
    end;
%    
    defoci(miIndex)=mi.ctf.defocus;
    relTimes(miIndex)=(mov.timeList(miIndex)-min(mov.timeList))*24;
    jpegName=[mi.jpegPath mi.baseFilename '.jpg'];
    %     print(jpegName,'-djpeg');
    epaScl=1/max(epaVals.epaRawMinusBg);
    % put this onto the drift plot
    if graphicsOn    %%%%
        
        mysubplot(2,4,2);
        imags(ctfImage);
        axis off;
        title(['\delta = ' num2str(mi.ctf.defocus,3) '\mum  astig = ' num2str(abs(mi.ctf.deltadef),3) '\mum']);
        
        mysubplot(4,2,5);
        plot(1./epaVals.resolution,[epaVals.ctfSim.^2 epaScl*epaVals.epaRawMinusBg 0*epaVals.ccc epaVals.ccc]);
        title(['Res limit: ' num2str(ctfVals.RES_LIMIT,3) 'Å']);
        
        mysubplot(4,2,7);
        plot(relTimes,defoci,'o');
        ylabel('Defocus, \mum');
        xlabel('Time, hr');
        drawnow;
        print(jpegName,'-djpeg');
        
    else
        %%%%%
        disDat.frameShifts.title=['Res limit: ' num2str(ctfVals.RES_LIMIT,3) 'Å'];
        s=struct;
        s.exec='mysubplot(242)';
        s.image=ctfImage;
        s.axis='off';
        s.title=['\delta = ' num2str(mi.ctf.defocus,3) '\mum ' ...
            ' astig = ' num2str(abs(mi.ctf.deltadef),3) '\mum'];
        disDat.ctfImage=s;
        %%%%%
        s=struct;
        s.exec='mysubplot(425)';
        s.xplot=1./epaVals.resolution;
        s.yplot=[epaVals.ctfSim.^2 epaScl*epaVals.epaRawMinusBg 0*epaVals.ccc epaVals.ccc];
        disDat.epa=s;
        
        %%%%%
        s=struct;
        s.exec='subplot(427)';
        s.xplot=relTimes;
        s.yplot=defoci;
        s.ylabel='Defocus \mum';
        s.xlabel='Time, hr';
        
        disDat.defocusPlot=s;
        
        %%%%
        disDat.print.eval=['print(' jpegName '''-djpeg'')'];
        save([mi.jpegPath mi.baseFilename 'dis.dat']);
        
    end;
    
    %     if graphicsOn
    %     DrawFigureFromData(disDat);
    %     drawnow;
    %     end;
    miName=[infoDir mi.baseFilename 'mi.txt'];
    WriteMiFile(mi,miName);
    disp(['Wrote ' miName]);
    disp(' ');
end;  % big while


function s=rtFindNextMovie(s,mode)
% modes: 1 latest; 2 scan from beginning; 3 from beginning, one pass.
% s is a struct containing nameList, timeList, where the last entry is the
% current file and its creation time.
% needs:
% s.doInit
% s.inPath
% s.foundAFile
% s.nameList --last entry is newest name
% s.timeList --last entry is newest time
% s.freeList --files that have not been operated on
% s.fileIndex -- index into the directory

s.foundAFile=false;

if numel(s.nameList)<1
    s=GetNewFiles(s);
end;

switch mode
    case 1  % just pick the last
        if ~any(s.freeList)
            s=GetNewFiles(s);
        end;
        if any(s.freeList)
            s.fileIndex=find(s.freeList,1,'last');
            s.foundAFile=true;
        end;
        
    case {2 3}  % starting from the beginning, work forward.
        if any(s.freeList)
            s.fileIndex=find(s.freeList,1,'first');
            s.foundAFile=true;
        elseif mode==3 % repeat
            s=GetNewFiles(s);
            if any(s.freeList)
                s.fileIndex=find(s.freeList,1,'first');
                s.foundAFile=true;
            end;
        end;
    otherwise
        error(['Unrecognized mode ' num2str(mode)]);
end;

if s.foundAFile
    s.freeList(s.fileIndex)=false;
end;

function s=GetNewFiles(s);
% to get initial set of files, call with s.nameList={}, s.timeList=[],
% s.fileIndex=0;
% to get additional files, call with s.fileIndex=numel(s.nameList)
numFiles=numel(s.nameList);
d=dir(s.inPath);
for i=1:numel(d)
    nameOk=numel(strfind(d(i).name,s.pattern));
    if ~d(i).isdir && nameOk
        if ~any(strcmp(d(i).name,s.nameList)) % a new name
            numFiles=numFiles+1;
            s.nameList{numFiles}=d(i).name;
            s.timeList(numFiles)=d(i).datenum;
            s.freeList(numFiles)=true;
        end;
    end;
end;
nNew=numFiles-numel(s.nameList);
disp([num2str(nNew) ' new file' ess(nNew) ' found.']);
end
end