% MergeImages.m

mode=0;  % Interactive
% mode=22;  % Yale F20 carbon images

% Default parameters
readIfPresent=0;    % Read the existing mi file, else create new ones.
WriteInfo=1;        % Store the new mi file
makeDirectories=1;  % if directories don't exist, create them.
ignoreOldCTFs=1;    % If a CTF fit already exists, don't use it to set the defocus range.
fitAndAlign=1;      % Fit CTF and align images before merging.
mcDS=2;             % Composite image downsample ratio
useCircMask=1;      % Apply circular antialias filter in meCombineImages
maxNim=inf;         % number of image sets to merge
kV=200;
iCamera=1;          % Gatan camera
cpe=16;             % Gatan camera default

removeCrystal=0;
overridePixA=0;
overrideDose=0;
doRemoveGoldParticles=0;

startImage=1;
endImage=inf;
flipFirst=0;

ctfOptions=struct;
ctfOptions.lowB=1;  % lower B factors
ctfOptions.spFromWholeMicrograph=1;
ctfOptions.kV=kV;

% Start loading the micrograph structure
mi=meCreateMicrographInfoStruct13;
mi.imagePath='Micrograph/';
mi.infoPath='Info/';
mi.procPath='Merged/';
mi.tempPath='Temp/';
mi.kV=kV;
mi.camera=iCamera;
mi.cpe=cpe;


switch mode
    
    case 0
        
        % Load micrographs from the file selector

        [fname, pathName]=uigetfile({'*.mrc';'*.tif'},'Select micrograph files','multiselect','on');
        if isnumeric(pathName) % File selection was cancelled
            return
        end;
        if ~iscell(fname)
            fname={fname};
        end;
        [rootPath, dirMicrograph]=ParsePath(pathName);
        cd(rootPath);
        
        mi.basePath=rootPath;
        mi.imagePath=dirMicrograph;
        
        %     Examine the filenames
        [pa,nm,ex]=fileparts(fname{1});
        
        if strcmp(ex,'.mrc')  % presumably a K2 dataset
            disp([num2str(numel(fname)) ' files selected.  First is:']);
            disp(fname{1});
            disp('K2 dataset assumed.');
            mi.cpe=1;
            mi.camera=5;
            j=2;
            [pa,nm2]=fileparts(fname{j});
            while num2str(nm2(end))>num2str(nm(end))
                j=j+1;
                nm=nm2;
                if numel(fname)>=j
                    [pa,nm2]=fileparts(fname{j});
                end;
            end;
            ni=j-1;
            disp([num2str(ni) ' exposures per set.']);
            rexp='.+\.mrc';              % e.g. xxx.tif
            infos={};
            mi0=mi;  % default info structure
            for i=1:numel(fname)/ni
                j=(i-1)*ni+1;
                mi=mi0;  % copy the generic struct
                [pa mi.baseFilename ext]=fileparts(fname{j});
                mi.imageFilenames{1}=fname{j};
                for k=2:ni
                    mi.imageFilenames{k}=fname{j-1+k};
                end;
                infoFilename=[mi.basePath mi.infoPath mi.baseFilename 'mi.mat'];
                if exist(infoFilename,'file') && readIfPresent % the file exists already
                    disp([infoFilename ' -loaded']);
                    load(infoFilename);  % load the existing one.
                end;
                infos{i}=mi;
            end;
            
        elseif strcmp(ex,'.tif')  % .tif files from CCD, presumably
            ni=3;
            ni=MyInput('Number of exposures per set? ',ni);
            mi.camera=MyInput('Camera (1=F20CCD) ? ',1);
            mi.cpe=MyInput('Counts per electron (SerialEM CCD=9)',9);
            overridePixA=MyInput('pixA?', 1.7);
            
            rexp='.+\.tif';              % e.g. xxx.tif
            %           rexp='.+_1_\d+\.mrc';  e.g. xxx_1_123.mrc
            %           mi.cpe=44;     % SerialEM
            infos=meScanHidekiFiles2(mi,3,readIfPresent,rexp);
            
        end;
        
        switch ni
            case 2
                initialDefoci=[1.5 10];
            case 3
                initialDefoci=[1.5 8 16];
            otherwise
                error('initialDefoci not set up');
        end;
        startImage=1;
        endImage=inf;
    
%         -------special cases -------
        
    case 15  % DE12 DSR data
        numExposures=3;
        initialDefoci=[2 7 15]; %
        mi.basePath='/Volumes/cryoEMdata/supervisor/data/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/';
        mi.imagePath='Images/';
        aliDir='../Aligned/';
        flipFirst=-1  % special case to match drift-tracking
        mi.procPath='Merge/';
        mi.infoPath='Info/';
        
        infos=meScanDEFiles2(mi,numExposures);
        iCamera=3;  % DE12
        ctfOptions.lowB=1;  % lower B factors
        ctfOptions.spFromWholeMicrograph=1;
        mi.cpe=70;  %%counts per electron??
        overridePixA=2.9
        overrideDose=10
        %              startImage=11
        
    case 20  % Brandeis Falcon
        %             Sequential files, triples, Brandeis F20
        iCamera=4;  % Brandeis Falcon
        mi.camera=iCamera;
        mi.basePath='/Volumes/TetraData/EMWork/Hideki/121108_Brandeis/1/';
        readIfPresent=0;
        
        rexp='.+_\d+\.tif';
        mi.cpe=88;    % TIFF files from Brandeis
        overridePixA=2.04;
        %           rexp='.+_1_\d+\.mrc';
        %           cpe=44;     % SerialEM
        infos=meScanHidekiFiles2(mi,3,readIfPresent,rexp);
        initialDefoci=[2 5 10];
        
    case 21 %             Sequential files, triples, Yale F20
        %         mi.basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/121118/PC1PC2liposome_SUP_14k/';
        %         mi.basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/121118/PC1PC2liposome_SUP_14k/';
        mi.basePath='/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/';
        %         mi.basePath='/Volumes/TetraData/EMWork/Hideki/130604vesicles/test@80k/';
        %         mi.basePath='/Volumes/TetraData/EMWork/Hideki/130604vesicles/PCPSChol721/';
        %         mi.basePath='/Volumes/TetraData/EMWork/Hideki/DNA-LND/data/FigC/';
        %         doRemoveGoldParticles=1
        mi.basePath='/EMWork/Hideki/130820/Glur2GFP_AuBSAx2dilution/';
        mi.imagePath='Micrograph/';
        %         overridePixA=1.7;
        %          overridePixA=1.1;
        rexp='.+_\d+\.tif';              % e.g. xxx_123.tif
        rexp='.+\.tif';              % e.g. xxx.tif
        %           rexp='.+_1_\d+\.mrc';  e.g. xxx_1_123.mrc
        %           mi.cpe=44;     % SerialEM
        infos=meScanHidekiFiles2(mi,3,readIfPresent,rexp);
        initialDefoci=[1.5 5 10];
        %         initialDefoci=[4 4 4];
        startImage=3;
        endImage=6;
        overridePixA=2;
        mi.cpe=9;
    case 21.1  % binned TIFF images.
        mcDS=1;
        mi.basePath='/EMWork/Hideki/130422part/LiposomeMerge2kx2k/';
        mi.imagePath='Micrograph/';
        overridePixA=3.4;
        rexp='.+\.tif';              % e.g. xxx.tif
        infos=meScanHidekiFiles2(mi,3,readIfPresent,rexp);
        initialDefoci=[1.5 5 10];
        startImage=1;
        endImage=inf;
        
    case 22 %             Sequential files, six, Yale F20
        mi.basePath='/Volumes/TetraData/EMWork/Hideki/130629carbon/';
        mi.imagePath='Micrograph/';
        %         overridePixA=1.7;
        %         rexp='.+_\d+\.tif';              % e.g. xxx_123.tif
        rexp='.+\.mrc';              % e.g. xxx.tif
        %           rexp='.+_1_\d+\.mrc';  e.g. xxx_1_123.mrc
        %           mi.cpe=44;     % SerialEM
        infos=meScanHidekiFiles2(mi,6,readIfPresent,rexp);
        initialDefoci=[.5 1 2 4 8 16];
        startImage=1;
        endImage=inf;
    case 23 %             Sequential files, from K2
        % Set all the mi values before calling ScanFiles
        %         mi.basePath='/Volumes/cryoEMdata/Hideki/SerialEM/F20/SPA/131005/AMPAR_liposome_x3_MPQX/';
        mi.basePath='/Users/fred/EMWork/Hideki/140105/Box_AMPAR_Liposome_26_10mMGlu_slot4/';
        % mi.basePath='/Users/fred/EMWork/Hideki/130929/Box_antagonist3_slot1/';
        
        mi.imagePath='Micrograph/';
        
        mi.cpe=1;
        mi.camera=5;
        rexp='.+\.mrc';              % e.g. xxx.tif
        %           rexp='.+_1_\d+\.mrc';  e.g. xxx_1_123.mrc
        infos=meScanHidekiFiles2(mi,2,readIfPresent,rexp);
        initialDefoci=[1.5 10];
        initialDefoci=[5 10];
        %     Re-do the images with bars
        startImage=1;
        endImage=inf;
        
        
end;  % switch

%% Main program %%

jpegPath=[mi.procPath 'jpeg/'];
filtPath=[mi.procPath 'filtered/'];
filtJPath=[filtPath 'jpeg/'];

modelSpectrum=CCDModelSpectrum2D(mi.camera);  % handle DE-12 or CCD
[nCx, nCy]=size(modelSpectrum);  % Use the non-normalized spectrum, as
% mePreWhiten normalizes it.

nmi=numel(infos);
figure(10);
SetGrayscale;
figure(1);
SetGrayscale;
%%
% Get all the directories we'll use for writing results
writePaths={[mi.basePath mi.infoPath]
    [mi.basePath mi.procPath]
    [mi.basePath jpegPath]
    [mi.basePath filtPath]
    [mi.basePath filtJPath]};
%if needed, create directories.  The imagePath must already exist.
for j=1:numel(writePaths);
    if ~DirectoryExists(writePaths{j})
        if makeDirectories
            mkdir(writePaths{j});
        else
            warning(['the path doesn''t exist: ' writePaths{j}]);
        end;
    end;
end;

tic
lastImage=min(nmi,endImage);
disp(['Merging image sets ' num2str(startImage) ' to ' num2str(lastImage)]);
for i=startImage:lastImage
    
    disp(['>>Image index ' num2str(i)]);
    mi=infos{i};  % pick up the structure
    
    %%
    % Read the images
    [m, mi.pixA, mi.doses]=meReadImagesNorm(mi,mi.cpe,flipFirst,overridePixA);
    nim=size(m,3);
    mi.weights=single(ones(1,nim));
    if overrideDose  % not valid from DE-12 images
        mi.doses=mi.doses*0+overrideDose;
    else
        doses=mi.doses
    end;
    pixA=mi.pixA
    sz=size(m);
    mi.imageSize=sz(1:2);
    % operate with the pre-whitening filter
    %     disp('mePreWhiten');
    m=mePreWhiten(m,modelSpectrum);
    %%
    if fitAndAlign
        % Fit the ctfs
        %         ncts=numel(mi.ctf);
        ncts=nim;
        if ~ignoreOldCTFs && ncts>0  % some CTF information present: read the defocus values
            for j=1:ncts
                initialDefoci(j)=mi.ctf(j).defocus;
            end;
        end;
        rawCTFs=struct([]);
        
        if doRemoveGoldParticles
            m2=RemoveGoldParticles(m);
        else
            m2=m;
        end;
        
        for j=1:nim
            [rawCTF]=meFitCTF(m2(:,:,j),mi.pixA,...
                initialDefoci(j),removeCrystal,ctfOptions);  % Actual CTF fitting
            if j==1
                rawCTFs=rawCTF;
            else
                rawCTFs(j)=rawCTF;
            end;
            defocus=[rawCTFs(j).defocus rawCTFs(j).deltadef]
            jctfName=[mi.basePath filtJPath mi.baseFilename num2str(j) '-spect.jpg'];
            print('-djpeg','-r150',jctfName);  % save the CTF window.
        end;
        
        %%
        % Align the images
        % use downsampling by 8!
        [mi.mergeMatrix, mi.ctf, Ps]=meAlignMultiExposures(m,mi.pixA,rawCTFs,mi.doses,8);
        [mi.mergeMatrix, mi.ctf, Ps]=meAlignMultiExposures(m,mi.pixA,rawCTFs,mi.doses,4,Ps);
        %             Ps
        %%
        % Store the CTF and alignment information
        if WriteInfo
            save([mi.basePath mi.infoPath mi.baseFilename 'mi.mat'],'mi');
        end;
        
    end;
    %%
    % Create the composite image mc
    [effctf, mc, mts]=meCombineImages(m,mi,mcDS,useCircMask);
    % Note that this effctf doesn't include ccd transfer function, but
    % we're not using it anyway.
    mc=-mc;  % go back to reversed contrast
    %             Make a mask indicating the valid edge of the image.
    mskSize=round(mi.imageSize(1)/4);
    msk=meMakeMergedImageMask(mskSize,mi.mergeMatrix);
    mi=meInsertMask(msk,mi,1);  % merge mask is the first one in the mask stack.
    %%
    % Write out the combined image and filtered version.
    procname=[mi.basePath mi.procPath mi.baseFilename];  % Use the processed image directory
    jpegname=[mi.basePath jpegPath mi.baseFilename];
    filtname=[mi.basePath filtPath mi.baseFilename];
    filtjname=[mi.basePath filtJPath mi.baseFilename];
    
    % For jpeg files we trim .1% of values from the intensity
    % histograms.
    WriteMRC(mc,mcDS*mi.pixA,[procname 'm.mrc']);
    WriteJpeg(mc,[jpegname 'm.jpg']);
    imacs(mc);
    disp(['Wrote merged images: ' procname 'm.mrc']);
    
    title(mi.baseFilename,'interpreter','none');
    
    % Show the image and spectrum
    set(0,'CurrentFigure',10);  % i.e. figure(10)
    QuickLookSpectrum(mc,mcDS*mi.pixA,[mi.baseFilename 'm.mrc']);
    
    %%
    if removeCrystal
        mcs=meRemoveCrystal(mc,mcDS*mi.pixA,6,1);  % crystal-subtracted image
        WriteMRC(mcs,mcDS*mi.pixA,[procname 'mc.mrc']);
        WriteJpeg(mcs,[procname 'mc.jpg']);
    end;
    
    %%
    if WriteInfo
        infoFullName=meSaveMiFile(mi);
        disp(['Info file written: ' infoFullName]);
    end;
    disp(' ');
    infos{i}=mi;  % put the info back into the array
    
end;
toc

disp([num2str(lastImage-startImage+1) ' image sets processed.']);
disp(' ');
