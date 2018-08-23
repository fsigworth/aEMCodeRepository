% rlMiFromStar.m
% From a micrographs_ctf.star file, read the image filenames and ctf parameters
% and create a set of mi files and padded, normalized micrograph files.
% We assume we're already in the relion project directory.  We create two
% folders, Info/ and Merged/ to contain the mi and micrograph files
% respectively.

basePath=pwd; % assume we're in the relion project directory

cameraIndex=5; % K2
cpe=16;  % counts per electron, assumed.
dose=40;
ds=4;  % downsampling factor for 'small' images

disp('Getting a star file');
[starName,starPath]=uigetfile('*.star');
if isnumeric(starPath) % user clicked Cancel
    return
end;
[names,dat]=ReadStarFile([starPath starName]);
d=dat{1}; % read only the first data block
nim=numel(d.rlnMicrographName);

basePath=AddSlash(basePath);
CheckAndMakeDir([basePath mi0.procPath]);
CheckAndMakeDir([basePath mi0.infoPath]);

mi0=meCreateMicrographInfoStruct14;
for i=1:nim
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
    else
        mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;
    
    
    %     Read the micrograph; pad and rescale it.
    mrcFilename=[mi.imagePath mi.imageFilenames{1}];
    if exist(mrcFilename,'file')
        [m,s]=ReadMRC(mrcFilename);
        n0=size(m);
        n=NextNiceNumber(n0,5,8);  % new size, e.g. 3840 x 3840
        mi.imageSize=n;
        
        %     Pad and scale the image
        me=mean(m(:));
        m1=Crop(m,n,0,me);
        m2=(m1-me)/(5*std(m1(:))); % arbitrary simple scaling.
        
        %     Write out processed images into the Merged folder.
        WriteMRC(m2,mi.pixA,[mi.procPath mi.baseFilename 'm.mrc']);
        WriteMRC(Downsample(m2,n/ds),mi.pixA*ds,[mi.procPath mi.baseFilename 'ms.mrc']);
        
        WriteMiFile(mi);
    else
        disp(['Micrograph file not found: ' mrcFilename]);
    end;
end;

