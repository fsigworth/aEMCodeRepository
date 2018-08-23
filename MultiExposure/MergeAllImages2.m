% MergeAllImages21.m

mode=20;  % Brandeis F20
mode=21;  % Yale F20

% Default parameters
fitAndAlign=1;      % Fit CTF and align.  0 means merge only, assuming that
WriteInfo=1;        % Should we store the new mi data
ignoreOldCTFs=1;    % If a CTF fit already exists, don't use it to set the defocus range.
mcDS=2;             % Composite image downsample ratio
useCircMask=0;      % Apply circular antialias filter in meCombineImages
maxNim=inf;         % number of image sets to merge
uiPath='';
iCamera=1;          % Gatan camera
kV=200;
removeCrystal=0;
makeDirectories=1;  % if directories don't exist, create them.
makeSingleImages=0;      % =1, 2 or 3 if that single image is to be outputted.
overridePixA=0;
overrideDose=0;
doFilter=0;
doRemoveGoldParticles=0;
% load or, if necessary, create the micrograph info structures

startImage=1;
endImage=inf;
% endImage=2;

cpe=17.5;  % Counts per electron
aliDir='';
flipFirst=0;

jpegPath='';
filtPath='';  % additional paths

while 1  % repeating for the case of interactive merging.
    
    switch mode
           
            
        case 15  % DE12 DSR data
            mi=meCreateMicrographInfoStruct11;
            numExposures=3;
            initialDefoci=[2 7 15]; %
            mi.basePath='/Volumes/cryoEMdata/supervisor/data/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/';
            
            mi.imagePath='Images/';
            aliDir='../Aligned/';
            flipFirst=-1  % special case to match drift-tracking
            mi.procPath='Merge/';
            mi.infoPath='Info/';
            
            basePath=mi.basePath;
            imagePath=mi.imagePath;
            infoPath=mi.infoPath;
            procPath=mi.procPath;
            jpegPath=[mi.procPath 'jpeg/'];
            filtPath=[mi.procPath 'filtered/'];
            filtJPath=[filtPath 'jpeg/'];
            
            infos=meScanDEFiles2(mi,numExposures);
            iCamera=3;  % DE12
            ctfOptions.lowB=1;  % lower B factors
            ctfOptions.spFromWholeMicrograph=1;
            cpe=70;  %%counts per electron??
            overridePixA=2.9
            overrideDose=10
            %              startImage=11
            mcDS=1  % no downsampling.
            
            
        case 20
            %             Sequential files, triples, Brandeis F20
            mi=meCreateMicrographInfoStruct11;
            mi.kV=200;
            iCamera=4;  % Brandeis Falcon
            mi.camera=iCamera;
            ctfOptions.kV=200;
            ctfOptions.spFromWholeMicrograph=1;
            mi.basePath='/Volumes/raid3/Hideki/Current_DATA/F20_Brandeis/121108/1/';
          mi.basePath='/Volumes/raid3/Hideki/Current_DATA/F20_Brandeis/121109/slot1/';
          mi.basePath='/Volumes/TetraData/EMWork/Hideki/121109_Brandeis/slot1/';
          mi.basePath='/Volumes/TetraData/EMWork/Hideki/121108_Brandeis/1/';
            mi.imagePath='Micrograph/';
%           mi.imagePath='sq06_1/';
            mi.infoPath='Info/';
            mi.procPath='Merged/';
            readIfPresent=0;
          
          rexp='.+_\d+\.tif';
          cpe=88;    % TIFF files from Brandeis
          overridePixA=2.04;

%           rexp='.+_1_\d+\.mrc';
%           cpe=44;     % SerialEM

            infos=meScanHidekiFiles2(mi,3,readIfPresent,rexp);
            initialDefoci=[2 5 10];
          
            basePath=mi.basePath;
            imagePath=mi.imagePath;
            infoPath=mi.infoPath;
            procPath=mi.procPath;
            jpegPath=[mi.procPath 'jpeg/'];
            filtPath=[mi.procPath 'filtered/'];
            filtJPath=[filtPath 'jpeg/'];

        case 21
            %             Sequential files, triples, Yale F20
            mi=meCreateMicrographInfoStruct11;
            mi.kV=200;
            iCamera=1;  % F20 US4000
            mi.camera=iCamera;
            ctfOptions.kV=200;
            ctfOptions.spFromWholeMicrograph=1;
          mi.basePath='/Volumes/TetraData/EMWork/Hideki/121110/slot2/';
          mi.basePath='/Volumes/TetraData/EMWork/Hideki/121110/slot3/';
          mi.basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/121118/PC1PC2liposome_SUP_14k/';
          mi.basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/121118/PC1PC2liposome_SUP_14k/';
            mi.basePath='/Volumes/TetraData/EMWork/Hideki/121119/10%chol_SUP_14k/';
            mi.basePath='/Volumes/TetraData/EMWork/Hideki/121121/AMPARliposome_apo_SUP_14kx6_Box_liposome_tmp3_slot1/';
mi.basePath='/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/';
mi.basePath='/Volumes/TetraData/EMWork/Hideki/121206p/AMPARliposome_antagonist_slot1/';

mi.basePath='/Volumes/TetraData/EMWork/Hideki/DNA-LND/data/FigC/';
          mi.basePath='/Users/fred/Desktop/Data/121210/10003/';



% cpe=9;  % serialEM ????????
makeSingleImages=0;
fitAndAlign=1;
startImage=1;
cpe=15;  % guess
overridePixA=1.7;
doRemoveGoldParticles=0;

% mi.basePath='/Users/fred/EMWorkCache/Hideki/121219p/';
% mi.basePath='/Volumes/TetraData/EMWork/Hideki/121225/BoxGII5_slot2/';  % bad merging example
%             ctfOptions.lowB=1;  % lower B factors
% startImage=8;
% endImage=8;

            mi.imagePath='Micrograph/';
%             mi.imagePath='Micrograph.sq02_1/';
            mi.infoPath='Info/';
            mi.procPath='Merged/';
            readIfPresent=0;
          
%             maxNim=2;
          
          rexp='.+_\d+\.tif';
%               overridePixA=1.7;
%             overrideDose=20;
            
%           rexp='.+_1_\d+\.mrc';
%            cpe=44;     % SerialEM

            infos=meScanHidekiFiles2(mi,3,readIfPresent,rexp);
            initialDefoci=[2 5 10];
            initialDefoci=[1.5 5 10];
          
            basePath=mi.basePath;
            imagePath=mi.imagePath;
            infoPath=mi.infoPath;
            procPath=mi.procPath;
            jpegPath=[mi.procPath 'jpeg/'];
            filtPath=[mi.procPath 'filtered/'];
            filtJPath=[filtPath 'jpeg/'];
            

            
        case 100 % interactive merging of a file pair or triple
            removeCrystal=0;
            findVesicles=0;
            [files uiPath]=uigetfile('*','Select an image set (low def comes first)','multiselect','on',uiPath);
            imagePath=uiPath;
            procPath=uiPath;
            infoPath=uiPath;
            mi=meCreateMicrographInfoStruct10;
            
            if isnumeric(files);
                return
            end;
            if ischar(files)
                mi.imageFilenames{1}=files;
                nim=1;
            else
                mi.imageFilenames=files;
                nim=numel(files);
            end;
            [pa1 nm ex]=fileparts(mi.imageFilenames{1});
            mi.baseFilename=nm;
            initialDefoci=[2 5 15];
            infos={};
            infos{1}=mi;
            startImage=1;
    end;  % switch
    
    %% Main program %%
    
    
    modelSpectrum=CCDModelSpectrum2D(iCamera);  % handle DE-12 or CCD
% We should normalize it; but right now we don't.  Maybe change this in version 122
%    modelSpectrum=modelSpectrum/max(modelSpectrum(:));  %%Normalize it.
ctfOptions=struct;
ctfOptions.kV=kV;  % 200 keV
    
    nmi=numel(infos);
    figure(10);
    SetGrayscale;
    figure(1);
    SetGrayscale;
    %%
    % Get all the directories we'll use for writing results
    writePaths={[basePath infoPath]
        [basePath procPath]
        [basePath jpegPath]
        [basePath filtPath]
        [basePath filtJPath]};
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
    disp(['Images ' num2str(startImage) ' to ' num2str(min(nmi,endImage))]);
    for i=startImage:min(nmi,endImage)
        
        disp(['>>Image index ' num2str(i)]);
        mi=infos{i};  % pick up the structure
        mi.basePath=basePath;
        mi.imagePath=imagePath;
        mi.procPath=procPath;
        mi.infoPath=infoPath;
        % Construct the complete image file names
        nim=numel(mi.imageFilenames);
        nim=min(maxNim,nim);
        fnames=cell(1,nim);
        set(0,'CurrentFigure',1);  % replacement for figure(1), but not superior...
        disp('Processing images:');
        for j=1:nim
            fnames{j}=[mi.basePath mi.imagePath mi.imageFilenames{j}];
            disp(mi.imageFilenames{j});
        end;
        %%
        % Read the images
        if numel(aliDir)>0
            name1=fnames{1};
            [pa nm ex]=fileparts(name1);
            fnames{1}=[AddSlash(pa) aliDir nm '-ali.mrc'];
            %             char(fnames)
        end;
        
        [m mi.pixA mi.doses]=meReadImages(fnames,cpe,flipFirst,overridePixA);
        
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
                print('-djpeg','-r0',jctfName);  % save the CTF window.
            end;
            
            %%
            % Align the images
            % use downsampling by 8!
            [mi.mergeMatrix mi.ctf Ps]=meAlignMultiExposures(m,mi.pixA,rawCTFs,mi.doses,8);
            [mi.mergeMatrix mi.ctf Ps]=meAlignMultiExposures(m,mi.pixA,rawCTFs,mi.doses,4,Ps);
%             Ps
            %%
            % Store the CTF and alignment information
            if WriteInfo
                save([basePath infoPath mi.baseFilename 'mi.mat'],'mi');
            end;
            
        end;
        %%
        % Create the composite image mc
            [effctf mc mts]=meCombineImages(m,mi.pixA,mi.ctf,mi.mergeMatrix,mi.doses,mcDS,useCircMask);
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
            infoFullName=[basePath infoPath mi.baseFilename 'mi.mat'];
            save(infoFullName,'mi');
            disp(['Info file written: ' infoFullName]);
        end;
        infos{i}=mi;  % put the info back into the array
        
    end;
    toc
    if mode ~=100
        return
    end;
end;  % while loop

disp([num2str(-startImage+min(nmi,endImage)) ' image sets processed.']);
disp(' ');
