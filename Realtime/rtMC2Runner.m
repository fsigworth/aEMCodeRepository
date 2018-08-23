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
doDoseWeighting=0;
countingImageSize=3840;
writeAlignedStack=1;
    patches=[3 3];
    sgpus='0 1 2 3';

    
    %% Start: read the first movie frame and get size, dose, perhaps pixA
%     if ~isa(mi.movieFilename,'cell')
%         mi.movieFilename={mi.movieFilename};
%     end;
    disp(' ');
    disp(['Processing ' mi.movieFilename{1}]);
    [m,s,ok]=ReadMovie([mi.moviePath mi.movieFilename{1}],1,1); % read one frame
    n0=size(m);
    ds=round(n0(1)/countingImageSize);  % get the downsampling factor
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
    if doDoseWeighting
        frameDose=mean(m(:))*ds^2/(mi.cpe*mi.pixA^2); % dose in e/A^2
        % we have to scale up because MC2 doesn't.
        mi.frameDose=frameDose*ones(nFrames,1);
        mi.doses=frameDose*nFrames;
        doseString=[' -FmDose ' num2str(mi.frameDose(1))];
    else
        doseString='';
    end;
    if writeAlignedStack
        stackString=' -OutStack 1 ';
    else
        stackString='';
    end;
    gainString='';
    mi.frameSets=[1 nFrames];
    
    sftBin=num2str(ds);
    nFrames=s.nz;
    CheckAndMakeDir(mi.tempPath,1);
    CheckAndMakeDir(mi.imagePath,1);
    
    movieName=mi.movieFilename{1};
    
    [pa,mvBaseName,ex]=fileparts(movieName);
    imageName=[mi.imagePath mi.imageFilenames{1}];
    
%     disp(['Reading ' movieName]);
%     mv=ReadMovie(movieName);
%     ex='.mrcs';
%     movieName=[mvBaseName '.mrcs'];
%     mi.movieFilename{1}=movieName;
%     disp(['Writing ' movieName]);
%     WriteMRC(mv,mi.pixA,[mvBaseName ex]);
    
    MC2Exec='MotionCor2_1.1.0-Cuda80';
    if strcmp(ex,'.tif')
        movieType='InTiff';
    else
        movieType='InMrc';
    end;
    
    % Create the execution script
    string=[MC2Exec ' -' movieType ' ' [mi.moviePath mi.movieFilename{1}]...
        ' -OutMRC ' imageName ...
        ' -LogFile ' [mi.tempPath mvBaseName] '-' ...
          gainString stackString...
        ' -Patch ' num2str(patches) ' -Gpu ' sgpus ' -Kv ' num2str(mi.kV) ...
        ' -FtBin ' sftBin ' -PixSize ' num2str(mi.pixA/ds) doseString ...
        ' >> ' [mi.tempPath 'MC2Out.txt'] ' 2>> ' [mi.tempPath 'MC2Out.err'] ];
    disp(string);

%     exf=fopen('temp/ExecMC2.sh','w');
%     fprintf(exf,'%s\n',string);
% fclose(exf);
% system('chmod a+x temp/ExecMC2.sh');
% system('temp/ExecMC2.sh');
% 
%     
    
    system(string);
    
    if doDoseWeighting
        % Correct the micrograph scaling
        dwName=[mi.tempPath mvBaseName '_DW.mrc']; % read the dose-weighted output
        [m,s]=ReadMRC(dwName);
        me=mean(m(:));
        m2=(m-me)*sqrt(nFrames);
        mOut=(Crop(m2,mi.imageSize)+me)*ds^2;  % pad the image, restore mean and scale up
        mi.imageFilenames={[mi.baseFilename 'ala.mrc']};
        outName=[mi.imagePath mi.imageFilenames{1}];
        WriteMRC(mOut,s.pixA,outName);  % write it back out.
        disp(['Corrected micrograph written: ' outName]);       
        % delete the original files
        system(['rm ' mi.tempPath mvBaseName '*.mrc']);
    end;
    %%
    % Read the alignment shifts
    if any(patches>1)
        suffix='-0-Patch-Full.log';
        headerLines=3;
    else
        suffix='-0-Full.log';
        headerLines=1;
    end;
    logName=[mi.tempPath mvBaseName suffix];
    if ~exist(logName,'file')
        disp(['Gctf log file not found: ' logName]);
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
