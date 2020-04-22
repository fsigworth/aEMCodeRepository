function rlStarToMiFiles(starName,pars)
% function mi=rlStarToMiFiles(starName,pars)
% function mi=rlStarToMiFiles(starCells,pars)
% From a micrographs_ctf.star file, read the image filenames and ctf parameters
% and create a set of mi files. We create the
% Info/ directory to contain the mi files.
% if pars. writImage is selected, Merged/ is created and contains padded and scaled
% micrograph files.
% If the output argument mi is given, only one mi file is generated, which will be the one from
% the line of index pars.startImage will be used.
% If the first argument starCells is a cell array, it contains the outputs from
% ReadStarFile() so the file doesn't have to be read again.
if nargin<1
    starName=''; % we may have to put up the file selector.
end;
if nargin<2
    pars=struct; % use all defaults.
end;
dpars=struct;

dpars.starName='';
dpars.basePath=pwd; % assume we're in the relion project directory
dpars.cameraIndex=5; % K2
dpars.cpe=0.8;  % counts per electron, for K2 counting mode.
% cameraIndex=6; % 'Falcon2' as we don't have info for Falcon 3 yet.
% cpe=64;
% processed by MotionCor2, this should be 0.2 I think.
dpars.dose=60; % Approx total movie dose in e/A^2
dpars.ds=8;  % downsampling factor for 'small' image
dpars.writeFullSize=0; % write out full-size *.m image
dpars.writeDownsampled=0;
dpars.writeMiFile=1;
dpars.startImage=1;
if nargout==0
    dpars.lastImage=1;
else
    dpars.lastImage=inf;
end;
dpars.BFactor=60;
dpars.imagePath=[];
dpars.newMicrographPath='';
dpars.readMicrographForScale=1;

pars=SetDefaultValues(dpars,pars,1); % 1 means check for undefined fieldnames.

cd(pars.basePath);

if isa(starName,'cell')
    names=starName{1};
    dat=starName{2};
else
    if numel(pars.starName)<1
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
disp([num2str(nLines) ' micrographs.']);

for i=1:nLines
    [ok,mi,m]=rlStarLineToMi(names,dat,i);
    if ok
        CheckAndMakeDir(mi.infoPath);
        miName=WriteMiFile(mi);
        disp([num2str(i) ' ' miName]);
        nsz=NextNiceNumber(mi.imageSize);
        imags(BinImage(Crop(m,nsz,0,mean(m(:))),4));
        title(mi.baseFilename,'interpreter','none');
        drawnow;
    end;
end;

% iOptics=0; % pointer to Optics block in star file
% iMicrographs=1;
% opt=struct;
% mic=struct; % This will be the main 'micrograph' struct.
% for i=1:numel(names)
%     nm=names{i};
%     if iOptics==0 && numel(strfind(nm,'optics'))>0
%         iOptics=i;
%     elseif iMicrographs==1 && numel(strfind(nm,'micrographs'))>0
%         iMicrographs=i;
%     end;
% end;
% if iOptics>0 % There is an optics block
%     opt=dat{iOptics};
% end;
% mic=dat{iMicrographs};
% 
% nim=numel(mic.rlnMicrographName);
% 
% % basePath=AddSlash(pars.basePath);
% 
% mi0=meCreateMicrographInfoStruct14;
% 
% if nargout==0  % we'll be writing an mi file
%     CheckAndMakeDir(mi0.infoPath,1);
% end;
% if pars.writeImages
%     CheckAndMakeDir(mi0.procPath);
% end;
% lastImage=min(nim,pars.lastImage);
% 
% disp(['Processing ' num2str(lastImage) ' micrographs']);
% for i=pars.startImage:lastImage
%     mi=mi0; % copy the default parameters
%     [imagePath,baseName,imageExtension]=fileparts(mic.rlnMicrographName{i});
%     if numel(pars.imagePath)>0 % We want to change this?
%         imagePath=pars.imagePath;
%     end;
%     if ~exist(imagePath,'dir')
%         error(['The image path ' imagePath ' wasn''t found']);
%     end;
%     imagePath=AddSlash(imagePath);
%     
%     mi.baseFilename=baseName;
%     mi.originalBasePath=AddSlash(pwd);
%     mi.basePath='';
%     mi.moviePath='';
%     mi.imageFilenames{1}=[baseName imageExtension];
%     mi.imagePath=imagePath;
%     
%     % pick up or calculate the pixel size
%     if isfield(mic,'rlnMicrographPixelSize') || isfield(opt,'rlnMicrographPixelSize')
%         mi.pixA=GetOptField('rlnMicrographPixelSize',i);
%     else
%         mi.pixA=1e4*mic.rlnDetectorPixelSize(i)/mic.rlnMagnification(i);
%     end;
% 
%     mi.doses=pars.dose;
%     mi.kV=GetOptField('rlnVoltage',i);
%     mi.camera=pars.cameraIndex;
%     mi.cpe=pars.cpe;
%     mi.weights=1;
%     
% %   Set up to use the original micrograph as the "processed" or "merged"
% %   image.
%     mi.procPath=imagePath;
% %     Compute the normalization we'll have to apply to make the AC image component
% %     approximately equal to fractional contrast. We'll estimate the variance
% %     by averaging the 1D power spectrum from 0.3 to 0.7 * Nyquist. Assuming
% %     an arbitrary image scaling by a, the est variance of a counting image should be
% %     a^2 * pixelDose, which will be (a * size of original pixel)^2 * doses(1),
% %     with a being the unknown scaling factor. (By
% %     the scaling used by MotionCor2, the "size of original pixel" should be a 
% %     superres pixel. The final fractional-contrast image will
% %     have the variance 1/pixelDose. We get it by scaling the raw image by
% %     1/sqrt(est variance).
% %     In the end we'll get the scaled micrograph as
% %     scaledImg = (m-mi.imageMedian)*mi.imageNormScale;
%     if pars.readMicrographForScale
%         m=ReadMRC([mi.imagePath mi.imageFilenames{1}]);
%         mi.imageSize=size(m);
%         sds=ceil(min(mi.imageSize)/256); % Downsampling factor for spectrum
%         mc=Crop(m,256*sds); % Grab a square region of the image.
%         sp1=RadialPowerSpectrum(mc,0,sds);
%         spN=numel(sp1);
%         spLims=round(spN*[.3 .7]);
%         estVar=sum(sp1(spLims(1):spLims(2)))/-diff(spLims);
%         mi.imageNormScale=1/sqrt(estVar);
%     else % we assume a true counting camera, binned from superres, and just estimate it.
%         mi.imageNormScale=2/(mi.pixA*sqrt(mi.doses(1))); % factor of 2 for 
% %            2-binned superres images.
%     end;
%     mi.imageMedian=median(m(:)); % best estimate we have of the main image mean.
%         
%     %     CTF parameters
%     mi.ctf=struct;
%     mi.ctf.defocus=(mic.rlnDefocusU(i)+mic.rlnDefocusV(i))/2e4;
%     mi.ctf.deltadef=(mic.rlnDefocusU(i)-mic.rlnDefocusV(i))/2e4;
%     mi.ctf.theta=mic.rlnDefocusAngle(i)*pi/180;
%     mi.ctf.Cs=GetOptField(rlnSphericalAberration,i);
%     mi.ctf.alpha=GetOptField(rlnAmplitudeContrast,i);
%     mi.ctf.B=pars.BFactor;
%     mi.ctf.lambda=EWavelength(mi.kV);
%     
%     mi.mergeMatrix=eye(3);
%     
%     % Use the Grant&Grigorieff damage model
%     if mi.kV>250
%         mi.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
%     else % 200 kV code
%         mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
%     end;
%     
% % %     %     Read the micrograph; pad and rescale it.
% % %     mrcFilename=[mi.imagePath mi.imageFilenames{1}];
% % %     if exist(mrcFilename,'file') && (writeFullSize || WriteDownsampled)
% % %         disp([num2str(i) ': Reading ' mrcFilename]);
% % %         [m,s]=ReadMRC(mrcFilename);
% % %         n0=size(m);
% % %         n=NextNiceNumber(n0,5,8);  % new size, e.g. 3840 x 3840
% % %         mi.imageSize=n;
% % %         
% % %         %     Pad and scale the image. We no longer calculate absolute
% % %         %     contrast, but simply scale according to the STD of the image.
% % %         me=mean(m(:));
% % %         m1=Crop(m,n,0,me);
% % %         m2=(m1-me)/(5*std(m1(:))); % arbitrary simple scaling, rather than absolute contrast.
% % %        
% % %         %     Write out processed images into the Merged folder.
% % %         if writeFullSize
% % %             disp(['Writing ' num2str(n) ' pixels: ' mi.procPath mi.baseFilename 'm.mrc']);
% % %             WriteMRC(m2,mi.pixA,[mi.procPath mi.baseFilename 'm.mrc']);
% % %         end;
% % %         if writeDownsampled
% % %             nd=n/ds;
% % %             m2d=Downsample(m2,nd);
% % %             disp(['Writing ' num2str(nd) ' pixels: ' mi.procPath mi.baseFilename 'ms.mrc']);
% % %             WriteMRC(m2d,mi.pixA*ds,[mi.procPath mi.baseFilename 'ms.mrc']);
% % %         end;
% % %     else
% % %         disp(['Micrograph file not written: ' mrcFilename]);
% % %         [~,s]=ReadMRC(mrcFilename,1,0); % Get just the header
% % %         n0=[s.nx s.ny];
% % %         n=NextNiceNumber(n0,5,8);  % new size, e.g. 3840 x 3840
% % %         mi.imageSize=n;
% % %     end;
% % %     if writeMiFile
% % %         disp(['Writing ' mi.baseFilename 'mi.txt']);
% % %         writeMiFile(mi);
% % %     end;
% end;
% 
% function val=GetOptField(fieldName,ind)
% % if mic.(fieldName) doesn't exist, pick it up from the opt structure.
% if isfield(mic,fieldName)
%     val=mic.(fieldName);
% else
%     iGrp=mic.rlnOpticsGroup(ind);
%     val=opt.(fieldName)(iGrp);
% end
% end
% 
% end
