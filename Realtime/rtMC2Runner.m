function mi=rtMC2Runner(mi,pars)
% ds is the downsampling factor, e.g. 2 for superres 8k stack to 4k output.
% The input mi file must have the following fields specified.
% basePath
% moviePath % typically should be absolute
% movieFilename
% gainRefName
% baseFilename
% imagePath
% imageFilenames
% imageSize
% tempPath % where we'll put log files
% pixA (replaced if zero)
% cpe
% kV

% For the returned miOut we fill in the additional fields
% imageSize
% frameDose
% frameSets
% doses
% damageModelCode
% frameShifts

% Example:
% mi.moviePath='movie_frames/';
% mi.movieFilename{1}='Aug10_14.58.56.tif';
% mi.gainRefName='CountRef.mrc';
% mi.baseFilename='Aug10_14.58.56';
% mi.imagePath='Micrograph/';
% mi.imageFilenames{1}=[mi.baseFilename 'c.mrc'];
% mi.tempPath='Log/';
% mi.pixA=0;
% mi.cpe=0.8;
% mi.kV=300;
dpars.throw=0;
dpars.trunc=0;
dpars.patches=[3 3];
dpars.writeAlignedStack=0;
dpars.defaultFrameDose=0;
if nargin<2
    pars=struct;
end;
pars=SetDefaultValues(dpars,pars);

deleteDWImage=1;
doDoseWeighting=1; % also computes frame doses
countingImageSize=3840;
    sgpus='0 1 2 3';

    
    %% Start: read the first movie frame and get size, dose, perhaps pixA
%     if ~isa(mi.movieFilename,'cell')
%         mi.movieFilename={mi.movieFilename};
%     end;
    disp(' ');
    disp(['Processing ' mi.movieFilename{1}]);
    [m,s,ok]=ReadMovie([mi.moviePath mi.movieFilename{1}],1,1); % read one frame
    if ~ok
        disp('Not a movie?');
        return
    end;
    n0=size(m);
    ds=max(1,round(n0(1)/countingImageSize));  % get the downsampling factor
    if ds>1
        disp(['Downsampling by ' num2str(ds)]);
    end;
    mi.imageSize=NextNiceNumber(n0/ds,5,8);
    nFrames=s.nz;
    if nFrames<2 % if it's not a movie, return mi=0.
        disp('Not a movie.');
        mi=0;
        return
    end;
    if mi.pixA==0
        mi.pixA=s.pixA*ds;
        disp(['Updated mi.pixA to ' num2str(mi.pixA)]);
    end;
    
    % Handle frame numbers and doses
    fr1=1+pars.throw;
    if fr1>nFrames
        pars.throw=0;
        disp('Throw was too large, set to zero.');
        fr1=1;
    end;
    
    frn=nFrames-pars.trunc; % Last frame no.
    nFramesUsed=frn-fr1+1;
    if nFramesUsed<1
        pars.trunc=0;
        frn=nFrames;
        disp('Trunc was too large, set to zero.');
    end;
    mi.frameSets=[fr1 frn];
    frameDose=mean(m(:))*ds^2/(mi.cpe*mi.pixA^2); % dose in e/A^2
    if pars.defaultFrameDose % if it's nonzero, use it.
        mi.cpe=mi.cpe*frameDose/pars.defaultFrameDose;
        frameDose=pars.defaultFrameDose;
    end;
    % we have to scale up because MC2 doesn't.
    mi.frameDose=frameDose*ones(nFrames,1);
    mi.doses=frameDose*nFramesUsed;
    
    if doDoseWeighting
        doseString=[' -FmDose ' num2str(mi.frameDose(1))];
    else
        doseString='';
    end;
    
    sftBin=num2str(ds);

    if pars.writeAlignedStack
        stackString=' -OutStack 1 ';
    else
        stackString='';
    end;
   
    if numel(mi.gainRefName)>1
       gainString=[' -Gain ' mi.gainRefName ' -RotGain ' num2str(mi.gainRefRot)];
   else
    gainString='';
   end;
   
    CheckAndMakeDir(mi.tempPath,1);
    CheckAndMakeDir(mi.imagePath,1);
    
    movieName=mi.movieFilename{1};
    
    [pa,mvBaseName,ex]=fileparts(movieName);
    tempImageName=[mi.tempPath mi.imageFilenames{1}];
    
%     disp(['Reading ' movieName]);
%     mv=ReadMovie(movieName);
%     ex='.mrcs';
%     movieName=[mvBaseName '.mrcs'];
%     mi.movieFilename{1}=movieName;
%     disp(['Writing ' movieName]);
%     WriteMRC(mv,mi.pixA,[mvBaseName ex]);
    
%     MC2Exec='MotionCor2_1.1.0-Cuda80';

MC2Exec='$RELION_MOTIONCOR2_EXECUTABLE';
%MC2Exec='MotionCor2.exe';
% MC2Exec='/ysm-gpfs/apps/software/MotionCor2/1.2.2-fosscuda-2018b/MotionCor2_1.2.2-Cuda92'
% MC2Exec='MotionCor2';
% MC2Exec='MotionCor2_1.2.2-Cuda92';
if strcmp(ex,'.tif')
        movieType='InTiff';
    else
        movieType='InMrc';
    end;
    
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
%     system(string);

    % Construct the execution script and run it.
    strings=cell(7,1);
    strings{1}=[MC2Exec ' -' movieType ' ' [mi.moviePath mi.movieFilename{1}] ' \'];
    strings{2}=['-OutMRC ' tempImageName ' \'];
    strings{3}=['-LogFile ' [mi.tempPath mvBaseName] '-' gainString stackString ' \'];
    strings{4}=['-Patch ' num2str(pars.patches) ' -Gpu ' sgpus ' -Kv ' num2str(mi.kV) ' \'];
    strings{5}=['-FtBin ' sftBin ' -PixSize ' num2str(mi.pixA/ds) doseString ' \'];
    strings{6}=['-Throw ' num2str(pars.throw) ' -Trunc ' num2str(pars.trunc) ' \'];
    strings{7}=['>> ' [mi.tempPath 'MC2Out.txt'] ' 2>> ' [mi.tempPath 'MC2Out.err'] ' \'];
    
    tempFile=[pars.tempDir 'MC2.sh'];
    exf=fopen(tempFile,'w');
    if exf<1
        error(['File could not be opened: ' tempFile]);
    end;
    for i=1:numel(strings)
        disp(strings{i});
        fprintf(exf,'%s\n',strings{i});
    end;
    disp(' ');
    fclose(exf);
    system(['chmod a+x ' tempFile]);
    system(tempFile);
   
    
    if ~exist(tempImageName,'file')
        disp('MC2: no output file.');
        mi=[];
        return
    end;
    
    %tempImageName=[mi.tempPath mi.imageFilenames{1}];
    %        mi.imageFilenames={[mi.baseFilename 'ala.mrc']};
    outName=[mi.imagePath mi.imageFilenames{1}]; % final micrograph name

    if doDoseWeighting
        % Correct the micrograph scaling
        [pa,nm,ext]=fileparts(tempImageName);
        dwName=[mi.tempPath nm '_DW.mrc']; % get the dose-weighted output
        disp(['Reading dose-weighted ' dwName]);
        [m,s]=ReadMRC(dwName);
        me=mean(m(:));
        m2=(m-me)*sqrt(nFrames); % scale up the AC part
        mOut=(Crop(m2,mi.imageSize)+me)*ds^2;  % pad the image, restore mean and scale up
        WriteMRC(mOut,s.pixA,outName);  % write it back out.
        disp(['Corrected micrograph written: ' outName]);       
        % delete the original files
        system(['rm ' tempImageName]);
        if deleteDWImage
            system(['rm ' dwName]);
        end;
    else
        disp(['No dose-weighting. Writing ' mi.imageFilenames{1}]);
        system(['mv ' tempImageName ' ' mi.imagePath mi.imageFilenames{1}]);
    end;
    %%
    % Read the alignment shifts
    if any(pars.patches>1)
        suffix='-0-Patch-Full.log';
        headerLines=3;
    else
        suffix='-0-Full.log';
        headerLines=1;
    end;
    logName=[mi.tempPath mvBaseName suffix];
    if ~exist(logName,'file')
        disp(['Gctf log file not found: ' logName]);
        disp('No shifts recorded.');
        mi.frameShifts{1}=zeros(mi.frameSets(2),2);
    else
        f=fopen(logName);
        for i=1:headerLines
            fgetl(f);
        end;
        vals=cell2mat(textscan(f,'%f%f%f'));
        fclose(f);
        mi.frameShifts{1}=vals(:,2:3);
    end;
    %     WriteMiFile(mi,allNames{fileIndex});
    
    %%
    
    %     plot(mi.frameShifts{1});
    %     grid on;
    %     title(mi.imageFilenames{1},'interpreter','none');
    %     drawnow;
    
end
