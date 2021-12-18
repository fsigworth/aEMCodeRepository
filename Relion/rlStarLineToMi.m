function [readOk,micFound,mi,nLines]=rlStarLineToMi(sNames,sDats,iLine,pars)
% function [ok, mi,img,nLines]=rlStarLineToMi(sNames,sDats,iLine,pars)
% Given [sNames,sDats]=ReadStarFile(name), create an mi struct from the
% line index iLine. Optionally returns the scaled micrograph img and the
% total number of lines in the star structures.
% Handles optics groups too.
% ok==false when there is an error in the star struct for that line.
% mi.imageSize==[0 0] in the case that micrograph couldn't be read.

defPars.changeImagePath=0;
defPars.checkImagePath=0;
defPars.writeMergedImage=0;
defPars.writeMergedSmall=0;
defPars.dose=50;
defPars.cameraIndex=5;
defPars.cpe=0.8;
defPars.BFactor=100;
defPars.noDamageModel=0;

if nargin<4
    pars=struct; % use all defaults.
end;

pars=SetDefaultValues(defPars,pars); % no checking for unused fields.

% Default outputs
readOk=false;
micFound=false;
mi=[];


% % Here are example par values (from rlStarToMiFiles).

% pars.basePath=pwd; % assume we're in the relion project directory
% dpars.cameraIndex=5; % K2
% dpars.cpe=0.8;  % counts per electron, 0.8 for K2 counting mode, but
% %  0.2 for superres image that is binned by MotionCor2.
% % ! For Falcon3: cameraIndex=6, I think cpe=16.
% dpars.dose=60; % Approx total movie dose in e/A^2. We have to guess this
% % because MotionCor2 scaling doesn't allow the total dose to be calculated.
% dpars.nFrames=40;
% dpars.BFactor=60; % Used by my CTF functions. Not critical.
% dpars.changeImagePath=0; % change from the path given in the star file
% dpars.imagePath='Micrographs/';
% dpars.defaultImageSize=[0 0]; % Value to insert if we can't read the
% % mrc file to get header information. If zeros, Vesicle_finding_GUI will
% % assign this upon reading the image.
% % dpars.readMicrographForScale=false; % If true, uses micrograph statistics
% % the actual image. Slower because each file is read.
% dpars.skipMissingMrcs=true; % Skip over any absent micrographs
% dpars.writeMiFile=1; % Write out each mi.txt file.
% dpars.writeMergedImage=1;
% dpars.writeMergedSmall=1;


% Pick up the parameters that might be special to an optics group
iOptics=0; % pointer to Optics block in star file
iMicrographs=0;
opt=struct;
% Find the first occurrences of '*optics' or '*micrographs'/'*particles' in sNames
for i=1:numel(sNames)
    nm=sNames{i};
    if iOptics==0 && numel(strfind(nm,'optics'))>0
        iOptics=i;
    elseif iMicrographs==0 && ...
        (numel(strfind(nm,'micrographs'))>0 || numel(strfind(nm,'particles'))>0)
        iMicrographs=i;
    end;
end;

iMicrographs=max(1,iMicrographs); % has to be at least 1

if iOptics>0 % There is an optics block
    opt=sDats{iOptics};
end;
mic=sDats{iMicrographs};

nLines=numel(mic.rlnMicrographName);
if iLine>nLines || iLine<1 % out of bounds, return nothing
    return
end;

mi0=meCreateMicrographInfoStruct14;

mi=mi0; % copy the default parameters
[imagePath,baseName,imageExtension]=fileparts(mic.rlnMicrographName{iLine});
if pars.changeImagePath
    imagePath=pars.imagePath;
end;

if pars.checkImagePath &&  numel(imagePath)>0 && ~exist(imagePath,'dir')
    error(['The image path ' imagePath ' wasn''t found']);
end;
imagePath=AddSlash(imagePath);

mi.baseFilename=baseName;
mi.originalBasePath=AddSlash(pwd);
mi.basePath='';
mi.moviePath='';
mi.imageFilenames{1}=[baseName imageExtension];
mi.imagePath=imagePath;

micName=[mi.imagePath mi.imageFilenames{1}];
micFound=exist(micName,'file');

if pars.writeMergedImage
    mi.procPath='Merged/';
end;
if pars.writeMergedSmall
    mi.procPath_sm='Merged_sm/';
end;

% pick up or calculate the pixel size
if isfield(mic,'rlnMicrographPixelSize') || isfield(opt,'rlnMicrographPixelSize') 
    mi.pixA=GetOptField('rlnMicrographPixelSize',iLine);
elseif isfield(opt,'rlnMicrographOriginalPixelSize')
    mi.pixA=GetOptField('rlnMicrographOriginalPixelSize',iLine);
else
    mi.pixA=1e4*mic.rlnDetectorPixelSize(iLine)/mic.rlnMagnification(iLine);
end;

mi.doses=pars.dose;
mi.kV=GetOptField('rlnVoltage',iLine);
mi.camera=pars.cameraIndex;
mi.cpe=pars.cpe;
mi.weights=1;

%     CTF parameters
mi.ctf=struct;
mi.ctf.defocus=(mic.rlnDefocusU(iLine)+mic.rlnDefocusV(iLine))/2e4;
mi.ctf.deltadef=(mic.rlnDefocusU(iLine)-mic.rlnDefocusV(iLine))/2e4;
mi.ctf.theta=mic.rlnDefocusAngle(iLine)*pi/180;
mi.ctf.Cs=GetOptField('rlnSphericalAberration',iLine);
mi.ctf.alpha=GetOptField('rlnAmplitudeContrast',iLine);
mi.ctf.B=pars.BFactor;
mi.ctf.lambda=EWavelength(mi.kV);

mi.mergeMatrix=eye(3);

if pars.noDamageModel
    mi.damageModelCode='1e6+f*0;'; % critical dose is huge constant
else
% Use the Grant&Grigorieff damage model
if mi.kV>250
    mi.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
else % 200 kV code
    mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
end;
end;

readOk=true; % We got the data from the star line.
return

    function val=GetOptField(fieldName,ind)
        % if mic.(fieldName) doesn't exist, pick it up from the opt structure.
        if isfield(mic,fieldName)
            val=mic.(fieldName)(ind);
        else
            iGrp=mic.rlnOpticsGroup(ind);
            val=opt.(fieldName)(iGrp);
        end
    end

end
