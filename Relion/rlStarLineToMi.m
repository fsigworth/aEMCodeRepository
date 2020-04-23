function [ok,mi,img,nLines]=rlStarLineToMi(sNames,sDats,iLine,pars)
% function [ok, mi,img,nLines]=rlStarLineToMi(sNames,sDats,iLine,pars)
% Given [smes,sDats]=ReadStarFile(name), create an mi struct from the
% line index iLine. Optionally returns the scaled micrograph img and the
% total number of lines in the star structures.
% Handles optics groups too.

if nargin<4
    pars=struct; % use all defaults.
end;

ok=0;
mi=[];
img=[];

% Here are example parameters.
% dpars.basePath=pwd; % assume we're in the relion project directory
% dpars.cameraIndex=5; % K2
% dpars.cpe=0.8;  % counts per electron, for K2 counting mode.
% % cameraIndex=6; % 'Falcon2' as we don't have info for Falcon 3 yet.
% % cpe=64;
% % processed by MotionCor2, this should be 0.2 I think.
% dpars.dose=60; % Approx total movie dose in e/A^2
% dpars.writeFullSize=0; % write out full-size *.m image
% dpars.writeDownsampled=0;
% dpars.ds=8;  % downsampling factor for 'small' image
% dpars.writeMiFile=0;
% dpars.BFactor=60;
% dpars.changeImagePath=1; % Where to find micrographs
% dpars.imagePath='';
% dpars.readMicrographForScale=1;
% pars=SetDefaultValues(dpars,pars,1); % 1 means check the fieldnames.

% Pick up the parameters that might be special to an optics group
iOptics=0; % pointer to Optics block in star file
iMicrographs=0;
opt=struct;
% Find the first occurrences of '*optics' or '*micrographs' in sNames
for i=1:numel(sNames)
    nm=sNames{i};
    if iOptics==0 && numel(strfind(nm,'optics'))>0
        iOptics=i;
    elseif iMicrographs==0 && numel(strfind(nm,'micrographs'))>0
        iMicrographs=i;
    end;
end;

if iMicrographs==0 % bad star structure
    warning('Bad star structure argument ''sDats''');
    return
end;

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
    if numel(imagePath)>0 && ~exist(imagePath,'dir')
        error(['The image path ' imagePath ' wasn''t found']);
    end;
    imagePath=AddSlash(imagePath);
    
    mi.baseFilename=baseName;
    mi.originalBasePath=AddSlash(pwd);
    mi.basePath='';
    mi.moviePath='';
    mi.imageFilenames{1}=[baseName imageExtension];
    mi.imagePath=imagePath;
    mi.procPath=imagePath; % look in the original images folder
    
    % pick up or calculate the pixel size
    if isfield(mic,'rlnMicrographPixelSize') || isfield(opt,'rlnMicrographPixelSize')
        mi.pixA=GetOptField('rlnMicrographPixelSize',iLine);
    else
        mi.pixA=1e4*mic.rlnDetectorPixelSize(iLine)/mic.rlnMagnification(iLine);
    end;

    mi.doses=pars.dose;
    mi.kV=GetOptField('rlnVoltage',iLine);
    mi.camera=pars.cameraIndex;
    mi.cpe=pars.cpe;
    mi.weights=1;

    if pars.readMicrographForScale || pars.skipIfNoImage
        micName=[mi.imagePath mi.imageFilenames{1}];
        if ~exist(micName,'file')
            return
        end;
    end;    
    
%   Set up to use the original micrograph as the "processed" or "merged"
%   image.
%     Compute the normalization we'll have to apply to make the AC image component
%     approximately equal to fractional contrast. We'll estimate the variance
%     by averaging the 1D power spectrum from 0.3 to 0.7 * Nyquist. Assuming
%     an arbitrary image scaling by a, the est variance of a counting image should be
%     a^2 * pixelDose, which will be (a * size of original pixel)^2 * doses(1),
%     with a being the unknown scaling factor. (By
%     the scaling used by MotionCor2, the "size of original pixel" should be a 
%     superres pixel. The final fractional-contrast image will
%     have the variance 1/pixelDose. We get it by scaling the raw image by
%     1/sqrt(est variance*pixelDose).
%     In the end we'll get the scaled micrograph as
%     scaledImg = (m-mi.imageMedian)*mi.imageNormScale;
    if pars.readMicrographForScale
        m=ReadMRC(micName);
        mi.imageSize=size(m);
        sds=floor(min(mi.imageSize)/256); % Downsampling factor for spectrum
        mc=single(Crop(m,256*sds)); % Grab a square region of the image.
        sp1=RadialPowerSpectrum(mc,0,sds);
        spN=numel(sp1);
        spLims=round(spN*[.3 .7]);
        estVar=sum(sp1(spLims(1):spLims(2)))/diff(spLims);
        mi.imageNormScale=1/(mi.pixA*sqrt(mi.doses*estVar));
        mi.imageMedian=median(m(:)); % best estimate we have of the main image mean.
        img=(single(m)-mi.imageMedian)*mi.imageNormScale;

    else % we assume a true counting camera, binned from superres, and 
%         just wing it for scale, don't compute median.
            mi.imageNormScale=1/(mi.cpe*mi.pixA^2*mi.doses(1));
    end;
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
    
    % Use the Grant&Grigorieff damage model
    if mi.kV>250
        mi.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else % 200 kV code
        mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;
    
    ok=1;
%     disp(iLine);
    return
% %     %     Read the micrograph; pad and rescale it.
% %     mrcFilename=[mi.imagePath mi.imageFilenames{1}];
% %     if exist(mrcFilename,'file') && (writeFullSize || WriteDownsampled)
% %         disp([num2str(iLine) ': Reading ' mrcFilename]);
% %         [m,s]=ReadMRC(mrcFilename);
% %         n0=size(m);
% %         n=NextNiceNumber(n0,5,8);  % new size, e.g. 3840 x 3840
% %         mi.imageSize=n;
% %         
% %         %     Pad and scale the image. We no longer calculate absolute
% %         %     contrast, but simply scale according to the STD of the image.
% %         me=mean(m(:));
% %         m1=Crop(m,n,0,me);
% %         m2=(m1-me)/(5*std(m1(:))); % arbitrary simple scaling, rather than absolute contrast.
% %        
% %         %     Write out processed images into the Merged folder.
% %         if writeFullSize
% %             disp(['Writing ' num2str(n) ' pixels: ' mi.procPath mi.baseFilename 'm.mrc']);
% %             WriteMRC(m2,mi.pixA,[mi.procPath mi.baseFilename 'm.mrc']);
% %         end;
% %         if writeDownsampled
% %             nd=n/ds;
% %             m2d=Downsample(m2,nd);
% %             disp(['Writing ' num2str(nd) ' pixels: ' mi.procPath mi.baseFilename 'ms.mrc']);
% %             WriteMRC(m2d,mi.pixA*ds,[mi.procPath mi.baseFilename 'ms.mrc']);
% %         end;
% %     else
% %         disp(['Micrograph file not written: ' mrcFilename]);
% %         [~,s]=ReadMRC(mrcFilename,1,0); % Get just the header
% %         n0=[s.nx s.ny];
% %         n=NextNiceNumber(n0,5,8);  % new size, e.g. 3840 x 3840
% %         mi.imageSize=n;
% %     end;
% %     if writeMiFile
% %         disp(['Writing ' mi.baseFilename 'mi.txt']);
% %         writeMiFile(mi);
% %     end;
% % end;

function val=GetOptField(fieldName,ind)
% if mic.(fieldName) doesn't exist, pick it up from the opt structure.
if isfield(mic,fieldName)
    val=mic.(fieldName);
else
    iGrp=mic.rlnOpticsGroup(ind);
    val=opt.(fieldName)(iGrp);
end
end

end
