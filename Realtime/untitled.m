% StarToMi.m
% Assuming we are in the project directory, we read the micrograph_ctf.star
% file and write mi.txt files into the new Info directory.  Also put
% normalized micrographs into the new Merged directory.

% Set these variables
microStarFile='all_micrographs_ctf.star';
cpe=20; % counts per electron
totalDose=60;  % e/A^2

% These are usually ok
mergeMode=3; % no phase flip; otherwise it would be 1.
bFactor=50;
infoDir='Info/';
procDir='Merged/';
jpegDir='Jpeg/';
ds=4; % downsampling factor for small image

CheckAndMakeDir(infoDir,1);
CheckAndMakeDir(procDir,1);
CheckAndMakeDir(jpegDir,1);

mi0=meCreateMicrographInfoStruct14;
mi0.procPath=procDir;
mi0.jpegPath=jpegDir;
mi0.moviePath='';
mi0.imagePath='';
mi0.doses=totalDose;
mi0.cpe=cpe;
mi0.camera=5;
mi0.weights=1;
mi0.mergeMode=mergeMode;
mi0.ctf=struct;
mi0.ctf.B=bFactor;

[names,dat]=ReadStarFile(microStarFile);
d=dat{1}; % get the first (only) set of data from the file.
figure(1);
axis equal;

for i=1:numel(d.rlnDefocusU) % Make an mi struct for each micrograph
    mi=mi0;
    [pa,mi.baseFilename]=fileparts(d.rlnMicrographName{i});
    mi.kV=d.rlnVoltage(i);
    if mi.kV>200
        mi.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else
        mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;
    
    mi.ctf.defocus=(d.rlnDefocusU(i) + d.rlnDefocusV(i))/2e4;
    mi.ctf.deltadef=(d.rlnDefocusU(i) - d.rlnDefocusV(i))/2e4;
    mi.ctf.theta=d.rlnDefocusAngle(i)*pi/180;
    mi.ctf.lambda=EWavelength(mi.kV);
    mi.ctf.Cs=d.rlnSphericalAberration(i);
    mi.ctf.alpha=d.rlnAmplitudeContrast(i);
    mi.ctf.ampFactor=1;
    mi.pixA=d.rlnDetectorPixelSize(i)/d.rlnMagnification(i)*1e4;

    img0=ReadMRC(d.rlnMicrographName{i});
    med=median(img0(:));
    n0=size(img0);
    n=NextNiceNumber(n0,4,8); % Force the image size to be a multiple of 8
    mi.imageSize=n;
    
    img=Crop(img0-med,n); % pad the image
    img=img/(totalDose*cpe*mi.pixA^2); % normalize to fractional contrast
    WriteMRC(img,mi.pixA,[procDir mi.baseFilename 'm.mrc']);
    imgSm=Downsample(img,n/ds);
    WriteMRC(imgSm,mi.pixA*ds,[procDir mi.baseFilename 'ms.mrc']);
    WriteJpeg(imgSm,[jpegDir mi.baseFilename 'ms.jpg']);
    imags(imgSm);
    title([num2str(i) '  ' procDir mi.baseFilename 'ms.mrc'],'interpreter','none');
    drawnow;
    infoName=WriteMiFile(mi);
    disp([infoName ' written.']);
end;
    
%     
%         %% Start: read the first movie frame and get size, dose, perhaps pixA
% %     if ~isa(mi.movieFilename,'cell')
% %         mi.movieFilename={mi.movieFilename};
% %     end;
%     disp(' ');
%     disp(['Processing ' mi.movieFilename{1}]);
%     [m,s,ok]=ReadMovie([mi.moviePath mi.movieFilename{1}],1,1); % read one frame
%     if ~ok
%         disp('Not a movie?');
%         return
%     end;
%     n0=size(m);
%     ds=round(n0(1)/countingImageSize);  % get the downsampling factor
%     if ds>1
%         disp(['Downsampling by ' num2str(ds)]);
%     end;
%     mi.imageSize=NextNiceNumber(n0/ds,5,8);
%     nFrames=s.nz;
%     if nFrames<2 % if it's not a movie, return mi=0.
%         disp('Not a movie.');
%         mi=0;
%         return
%     end;
%     if mi.pixA==0
%         mi.pixA=s.pixA*ds;
%         disp(['Updated mi.pixA to ' num2str(mi.pixA)]);
%     end;
%     
%     % Handle frame numbers and doses
%     fr1=1+pars.throw;
%     if fr1>nFrames
%         pars.throw=0;
%         disp('Throw was too large, set to zero.');
%         fr1=1;
%     end;
%     
%     frn=nFrames-pars.trunc; % Last frame no.
%     nFramesUsed=frn-fr1+1;
%     if nFramesUsed<1
%         pars.trunc=0;
%         frn=nFrames;
%         disp('Trunc was too large, set to zero.');
%     end;
%     mi.frameSets=[fr1 frn];
%     frameDose=mean(m(:))*ds^2/(mi.cpe*mi.pixA^2); % dose in e/A^2
%     % we have to scale up because MC2 doesn't.
%     mi.frameDose=frameDose*ones(nFrames,1);
%     mi.doses=frameDose*nFramesUsed;
%     
%     if doDoseWeighting
%         doseString=[' -FmDose ' num2str(mi.frameDose(1))];
%     else
%         doseString='';
%     end;
%     
%     sftBin=num2str(ds);
% 
%     if pars.writeAlignedStack
%         stackString=' -OutStack 1 ';
%     else
%         stackString='';
%     end;
%    
%     if numel(mi.gainRefName)>1
%        gainString=[' -Gain ' mi.gainRefName ' -RotGain ' num2str(mi.gainRefRot)];
%    else
%     gainString='';
%    end;
%    
%     CheckAndMakeDir(mi.tempPath,1);
%     CheckAndMakeDir(mi.imagePath,1);
%     
%     movieName=mi.movieFilename{1};
%     
%     [pa,mvBaseName,ex]=fileparts(movieName);
%     tempImageName=[mi.tempPath mi.imageFilenames{1}];
%     
% %     disp(['Reading ' movieName]);
% %     mv=ReadMovie(movieName);
% %     ex='.mrcs';
% %     movieName=[mvBaseName '.mrcs'];
% %     mi.movieFilename{1}=movieName;
% %     disp(['Writing ' movieName]);
% %     WriteMRC(mv,mi.pixA,[mvBaseName ex]);
%     
% %     MC2Exec='MotionCor2_1.1.0-Cuda80';
% 
% MC2Exec='$RELION_MOTIONCOR2_EXECUTABLE';
%     if strcmp(ex,'.tif')
%         movieType='InTiff';
%     else
%         movieType='InMrc';
%     end;
%     
%     % Create the execution script
%     string=[MC2Exec ' -' movieType ' ' [mi.moviePath mi.movieFilename{1}]...
%         ' -OutMRC ' tempImageName ...
%         ' -LogFile ' [mi.tempPath mvBaseName] '-' ...
%           gainString stackString ...
%         ' -Patch ' num2str(pars.patches) ' -Gpu ' sgpus ' -Kv ' num2str(mi.kV) ...
%         ' -FtBin ' sftBin ' -PixSize ' num2str(mi.pixA/ds) doseString ...
%         ' -Throw ' num2str(pars.throw) ' -Trunc ' num2str(pars.trunc) ...
%         ' >> ' [mi.tempPath 'MC2Out.txt'] ' 2>> ' [mi.tempPath 'MC2Out.err'] ];
%     disp(string);
% 
% %     exf=fopen('temp/ExecMC2.sh','w');
% %     fprintf(exf,'%s\n',string);
% % fclose(exf);
% % system('chmod a+x temp/ExecMC2.sh');
% % system('temp/ExecMC2.sh');
% % 
% %     
%     
%     system(string);
%     
%     if ~exist(tempImageName,'file')
%         disp('MC2: no output file.');
%         mi=[];
%         return
%     end;
%     
%     %tempImageName=[mi.tempPath mi.imageFilenames{1}];
%     %        mi.imageFilenames={[mi.baseFilename 'ala.mrc']};
%     outName=[mi.imagePath mi.imageFilenames{1}]; % final micrograph name
% 
%     if doDoseWeighting
%         % Correct the micrograph scaling
%         [pa,nm,ext]=fileparts(tempImageName);
%         dwName=[mi.tempPath nm '_DW.mrc']; % get the dose-weighted output
%         disp(['Reading dose-weighted ' dwName]);
%         [m,s]=ReadMRC(dwName);
%         me=mean(m(:));
%         m2=(m-me)*sqrt(nFrames); % scale up the AC part
%         mOut=(Crop(m2,mi.imageSize)+me)*ds^2;  % pad the image, restore mean and scale up
%         WriteMRC(mOut,s.pixA,outName);  % write it back out.
%         disp(['Corrected micrograph written: ' outName]);       
%         % delete the original files
%         system(['rm ' tempImageName]);
%         if deleteDWImage
%             system(['rm ' dwName]);
%         end;
%     else
%         disp(['No dose-weighting. Writing ' mi.imageFilenames{1}]);
%         system(['mv ' tempImageName ' ' mi.imagePath mi.imageFilenames{1}]);
%     end;
%     %%
%     % Read the alignment shifts
%     if any(pars.patches>1)
%         suffix='-0-Patch-Full.log';
%         headerLines=3;
%     else
%         suffix='-0-Full.log';
%         headerLines=1;
%     end;
%     logName=[mi.tempPath mvBaseName suffix];
%     if ~exist(logName,'file')
%         disp(['Gctf log file not found: ' logName]);
%         disp('No shifts recorded.');
%         mi.frameShifts{1}=zeros(mi.frameSets(2),2);
%     else
%         f=fopen(logName);
%         for i=1:headerLines
%             fgetl(f);
%         end;
%         vals=cell2mat(textscan(f,'%f%f%f'));
%         fclose(f);
%         mi.frameShifts{1}=vals(:,2:3);
%     end;
%     %     WriteMiFile(mi,allNames{fileIndex});
%     
%     %%
%     
%     %     plot(mi.frameShifts{1});
%     %     grid on;
%     %     title(mi.imageFilenames{1},'interpreter','none');
%     %     drawnow;
%     
% 
% 
%     
%     
%         m=ReadMRC([mi.imagePath mi.imageFilenames{1}]);
%     n0=size(m);
%     n=NextNiceNumber(n0,5,8);
%     m1=Crop(m,n,0,mean(m(:))); % expand to nice size
%     m1s=Downsample(m1,n/4);
% 
%     
%     
%     mi.ctf=struct;
% mi.ctf.defocus=(ctfVals.Defocus_U + ctfVals.Defocus_V)/2e4;
% mi.ctf.deltadef=(ctfVals.Defocus_U - ctfVals.Defocus_V)/2e4;
% mi.ctf.theta=ctfVals.Angle*pi/180;
% if pars.phasePlate
%     mi.ctf.phi=ctfVals.Phase_shift*pi/180;
% end;
% mi.ctf.alpha=pars.alpha;
% mi.ctf.lambda=EWavelength(mi.kV);
% mi.ctf.Cs=pars.Cs;
% mi.ctf.ampFactor=1;
% mi.ctf.ccc=ctfVals.CCC;
% mi.ctf.resLimit=ctfVals.RES_LIMIT;
% mi.ctf.estB=ctfVals.B_FACTOR;
% mi.ctf.B=pars.B;
% 
%     
    