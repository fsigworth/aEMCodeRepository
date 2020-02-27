% rlStarToMiFiles.m
% From a micrographs_ctf.star file, read the image filenames and ctf parameters
% and create a set of mi files and padded, normalized micrograph files.
% We assume we're already in the relion project directory.  We create two
% folders, Info/ and Merged/ to contain the mi and micrograph files
% respectively.

basePath=pwd; % assume we're in the relion project directory

cameraIndex=5; % K2
% cpe=0.8;  % counts per electron, for K2 counting mode. For super-res data
cameraIndex=6; % 'Falcon2' as we don't have info for Falcon 3 yet.
cpe=64;
% processed by MotionCor2, this should be 0.2 I think.
dose=60; % Approx total movie dose in e/A^2
ds=8;  % downsampling factor for 'small' images
writeFullSize=1; % write out full-size *.m image
writeDownsampled=1;
WriteMiFile=1;

disp('Getting a star file');
[starName,starPath]=uigetfile('*.star');
if isnumeric(starPath) % user clicked Cancel
    return
end;
%%
[names,dat]=ReadStarFile([starPath starName]);
d=dat{1}; % read only the first data block
nim=numel(d.rlnMicrographName);

basePath=AddSlash(basePath);

mi0=meCreateMicrographInfoStruct14;

CheckAndMakeDir([basePath mi0.procPath]);
CheckAndMakeDir([basePath mi0.infoPath]);
lastImage=nim;
% lastImage=1;
disp(['Processing ' num2str(lastImage) ' micrographs']);
for i=1:lastImage
    mi=mi0;
    [imagePath,baseName,ex]=fileparts(d.rlnMicrographName{i});
%     imagePath='';
    mi.baseFilename=baseName;
    mi.originalBasePath=AddSlash(basePath);
    mi.basePath=mi.originalBasePath;
    mi.moviePath='';
    mi.imagePath=AddSlash(imagePath);
    mi.imageFilenames{1}=[baseName ex];
    mi.pixA=1e4*d.rlnDetectorPixelSize(i)/d.rlnMagnification(i);
    mi.doses=dose;
    mi.kV=d.rlnVoltage(i);
    mi.camera=cameraIndex;
    mi.cpe=cpe;
    
    mi.weights=1;
    
    %     CTF parameters
    mi.ctf=struct;
    mi.ctf.defocus=(d.rlnDefocusU(i)+d.rlnDefocusV(i))/2e4;
    mi.ctf.deltadef=(d.rlnDefocusU(i)-d.rlnDefocusV(i))/2e4;
    mi.ctf.theta=d.rlnDefocusAngle(i)*pi/180;
    mi.ctf.Cs=d.rlnSphericalAberration(i);
    mi.ctf.alpha=d.rlnAmplitudeContrast(i);
    mi.ctf.B=40;
    mi.ctf.lambda=EWavelength(mi.kV);
    
    mi.mergeMatrix=eye(3);
    
    % Use the G&G damage model
    if mi.kV>250
        mi.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else % 200 kV code
        mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;
    
    
    %     Read the micrograph; pad and rescale it.
    mrcFilename=[mi.imagePath mi.imageFilenames{1}];
    if exist(mrcFilename,'file') && (writeFullSize || WriteDownsampled)
        disp([num2str(i) ': Reading ' mrcFilename]);
        [m,s]=ReadMRC(mrcFilename);
        n0=size(m);
        n=NextNiceNumber(n0,5,8);  % new size, e.g. 3840 x 3840
        mi.imageSize=n;
        
        %     Pad and scale the image. We no longer calculate absolute
        %     contrast, but simply scale according to the STD of the image.
        me=mean(m(:));
        m1=Crop(m,n,0,me);
        m2=(m1-me)/(5*std(m1(:))); % arbitrary simple scaling, rather than absolute contrast.
       
        %     Write out processed images into the Merged folder.
        if writeFullSize
            disp(['Writing ' num2str(n) ' pixels: ' mi.procPath mi.baseFilename 'm.mrc']);
            WriteMRC(m2,mi.pixA,[mi.procPath mi.baseFilename 'm.mrc']);
        end;
        if writeDownsampled
            nd=n/ds;
            m2d=Downsample(m2,nd);
            disp(['Writing ' num2str(nd) ' pixels: ' mi.procPath mi.baseFilename 'ms.mrc']);
            WriteMRC(m2d,mi.pixA*ds,[mi.procPath mi.baseFilename 'ms.mrc']);
        end;
    else
        disp(['Micrograph file not written: ' mrcFilename]);
        [~,s]=ReadMRC(mrcFilename,1,0); % Get just the header
        n0=[s.nx s.ny];
        n=NextNiceNumber(n0,5,8);  % new size, e.g. 3840 x 3840
        mi.imageSize=n;
    end;
    if WriteMiFile
        disp(['Writing ' mi.baseFilename 'mi.txt']);
        WriteMiFile(mi);
    end;
end;

