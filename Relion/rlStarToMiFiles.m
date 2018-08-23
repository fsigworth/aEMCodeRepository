% rlStarToMiFiles
% Takes the star file created by Gctf and makes a set of mi files.

% Currently set up to read Jesper's Piezo micrographs, e.g.
% Micrograph/piezo_560_2xsum_DW.mrc
% and the file Micrograph/micrographs_ctf.star

lp.cameraIndex=5; % K2
lp.cpe=0.2;  % should be 0.2 for superres movie binned by MC2
lp.bFactor=20;
lp.mergeMode=3;  % single exposure, no phase-flip
baseNameTrumcString='_DW';
imageNameAddString='_DW';
defaultAlpha=.02;

writeMiFiles=1;

rootPath=AddSlash(pwd);
micrographPath='Micrograph/';
moviePath='movie_frames/';
% refName='movie_frames/SuperRef_Aug10_14.58.56.dm4';
refName='';
dirInfo='Info/';
starName=[micrographPath 'micrographs_ctf.star'];

disp(['Reading ' starName]);
[blockNames,blockData,ok]=ReadStarFile(starName);
if ~ok
    error('Invalid star file.');
end;
%%
dat=blockData{1};
nmi=numel(dat.rlnMicrographName);

CheckAndMakeDir(dirInfo,1);

% Get the point at which we truncate the micrograph name to make the
% basename.  We examine the first micrograph name
p=strfind(dat.rlnMicrographName{1},'_DW'); % if _DW is present, truncate there.
if numel(p)<1                              % otherwise, just remove the extension.
    [pa,nm,ex]=fileparts(dat.rlnMicrographName{1});
    p=strfind(dat.rlnMicrographName{1},ex);
end;
baseNameEnd=p(1)-1;

% Create the basic mi structure
mi0=meCreateMicrographInfoStruct14;  % Get the mi file template
    % Put in the invariant parameters
mi0.camera=lp.cameraIndex;
mi0.cpe=lp.cpe;
mi0.originalBasePath=rootPath;
mi0.basePath=mi0.originalBasePath;
mi0.gainRefName=refName;
mi0.infoPath=dirInfo;
mi0.ctf=struct;
mi0.imagePath=micrographPath;
mi0.moviePath=moviePath;
mis=cell(nmi,1);
miNames=cell(nmi,1);

% Create mi files for the micrographs.
iStart=463;
iEnd=472;
for i=iStart:iEnd
%     prefix=sprintf('%05d_',i);
    prefix='';
    mi=mi0;
    
    name=dat.rlnMicrographName{i};
    [pa,nm,ex]=fileparts(name);
    imageName=[nm '_DW' ex];
    
%     truncFilename=dat.rlnMicrographName{i}(1:baseNameEnd);
%     truncFilename=name(1:baseNameEnd);
    truncFilename=nm;
    mi.baseFilename=[prefix truncFilename];
   % mi.imageFilenames=dat.rlnMicrographName(i);  % a cell scalar
    mi.imageFilenames={imageName};  % a cell scalar
    [~,s]=ReadMRC([mi.imagePath mi.imageFilenames{1}]);
    mi.imageSize=[s.nx s.ny];
    mi.kV=dat.rlnVoltage(i);
    mi.pixA=dat.rlnDetectorPixelSize(i)/dat.rlnMagnification(i)*1e4;
    mi.mergeMatrix=[1 0 0; 0 1 0; 0 0 0];
    defU=dat.rlnDefocusU(i);
    defV=dat.rlnDefocusV(i);
    mi.ctf.lambda=EWavelength(mi.kV);
    mi.ctf.defocus=(defU+defV)/2e4;
    mi.ctf.deltadef=(defU-defV)/2e4;
    mi.ctf.theta=dat.rlnDefocusAngle(i)*pi/180;
    if defaultAlpha
        mi.ctf.alpha=defaultAlpha;
    else
    mi.ctf.alpha=dat.rlnAmplitudeContrast(i);
    end;
    if isfield(dat,'rlnPhaseShift')
        mi.ctf.phi=dat.rlnPhaseShift(i)*pi/180;
    end;
%     mi.ctf.maxres=dat.rlnCtfMaxResolution(i);
    mi.ctf.Cs=dat.rlnSphericalAberration(i);
    mi.ctf.B=lp.bFactor;
    mi.ctf.ampFactor=1;
    mi.ctf.mergeMode=lp.mergeMode;
    if mi.kV>200
        mi.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else
        mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;
    miNames{i}=[mi.infoPath mi.baseFilename 'mi.txt'];
    if writeMiFiles
        WriteMiFile(mi,miNames{i});
    end;
    disp(miNames{i});
    mis{i}=mi;
end;
