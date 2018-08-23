% MergeAllImages21.m

mode=20;  % Brandeis F20

FitAndAlign=1;      % Fit CTF and align.  0 means merge only, assuming that
overwrite=1;        % That is, should we create new mi files.
WriteInfo=1;        % Should we store the new mi data
ignoreOldCTFs=1;    % If a CTF fit already exists, don't use it to set the defocus range.
findVesicles=0;
mcDS=2;             % Composite image downsample ratio
kWiener=.04;        % Wiener filter constant
fc=.1;              % Lowpass filter corner frequency applied to Wiener
useCircMask=0;      % Apply circular antialias filter in meCombineImages
maxNim=inf;
uiPath='';
iCamera=1;          % Gatan camera
removeCrystal=0;
ctfOptions=struct;
ctfOptions.kV=200;  % 200 keV
makeDirectories=1;  % if directories don't exist, create them.
singleImage=0;      % =1, 2 or 3 if that single image is to be outputted.
overridePixA=0;
overrideDose=0;
doFilter=0;
% load or, if necessary, create the micrograph info structures

startImage=1;
endImage=inf;

cpe=17.5;  % Counts per electron
aliDir='';
flipFirst=0;

jpegPath='';
filtPath='';  % additional paths

while 1  % repeating for the case of interactive merging.
    
    switch mode
           
            
        case 15  % DE12 DSR data
            mi=meCreateMicrographInfoStruct11;
            numExposures=3  % single exposure
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
%           mi.basePath='/Volumes/TetraData/EMWork/Hideki/121109_Brandeis/slot3/';
%           mi.basePath='/Volumes/TetraData/EMWork/Hideki/121108_Brandeis/1/';
%           mi.basePath='/Volumes/TetraData/EMWork/Hideki/121110/slot2/';
            mi.imagePath='Micrograph/';
           mi.imagePath='sq06_1/';
            mi.infoPath='Info/';
            mi.procPath='Merged/';
            readIfPresent=0;
          
          rexp='.+_\d+\.tif';
          cpe=88;    % TIFF files
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
    doFilter
    
    %% Main program %%
    
    
    modelSpectrum=CCDModelSpectrum2D(iCamera);  % handle DE-12 or CCD
    
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
        sz=size(m);
        mi.imageSize=sz(1:2);
        % operate with the pre-whitening filter
        %     disp('mePreWhiten');
        m=mePreWhiten(m,modelSpectrum);
        %%
        if FitAndAlign
            % Fit the ctfs
            %         ncts=numel(mi.ctf);
            ncts=nim;
            if ~ignoreOldCTFs && ncts>0  % some CTF information present: read the defocus values
                for j=1:ncts
                    initialDefoci(j)=mi.ctf(j).defocus;
                end;
            end;
            rawCTFs=struct([]);
            for j=1:nim
                [rawCTF]=meFitCTF(m(:,:,j),mi.pixA,...
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
            [mi.mergeMatrix mi.ctf]=meAlignMultiExposures(m,mi.pixA,rawCTFs,mi.doses,8);
            %%
            % Store the CTF and alignment information
            if WriteInfo
                save([basePath infoPath mi.baseFilename 'mi.mat'],'mi');
            end;
            
        end;
        %%
        % Create the composite image mc
        if singleImage
            for j=1:size(m,3)
            n=size(m);
            n(3)=[];
            mc0=AffineTransform(m(:,:,j),mi.mergeMatrix(:,:,j));
            mc1=Downsample(mc0,n/2);
            theCTF=CTF(n/2,mi.pixA*2,mi.ctf(j));
            mc=-real(ifftn(fftn(mc1).*ifftshift(sign(theCTF))));
            effctf=abs(theCTF);
            images(:,:,j)=mc;
            ctfs(:,:,j)=effctf;
            end;
            outName=[mi.procPath mi.baseFilename 'singles.mat'];
            save([mi.basePath outName],'images','ctfs');
            disp(['Wrote ' outName]);
        return
%             mcf=mc;
        else
            [effctf mc mts]=meCombineImages(m,mi.pixA,mi.ctf,mi.mergeMatrix,mi.doses,mcDS,useCircMask);
            mc=-mc;  % go back to reversed contrast
        end;
        %     [effctf mc]=rsCombineImages(m,mi,mcDS);
        %%
        % Filter the combined image
        if doFilter
            mcf=meFilterImage(mc,effctf,mcDS*mi.pixA,kWiener,fc);  % Filter the combined image
        end;
        %         lfAmp=0.2;
        %         mcf=meFilterImageFlatLF(mc,mi,lfAmp);
        
        % Write out the combined image and filtered version.
        procname=[mi.basePath mi.procPath mi.baseFilename];  % Use the processed image directory
        jpegname=[mi.basePath jpegPath mi.baseFilename];
        filtname=[mi.basePath filtPath mi.baseFilename];
        filtjname=[mi.basePath filtJPath mi.baseFilename];
   
        if singleImage  % make the filename ...img2m.mrc etc.
            procname=[procname 'img' num2str(singleImage)];
        end;
        % For jpeg files we trim .1% of values from the intensity
        % histograms.
        WriteMRC(mc,mcDS*mi.pixA,[procname 'm.mrc']);
        imwrite(uint8(imscale(rot90(mc),256,1e-3)),[jpegname 'm.jpg']);
        if doFilter
            WriteMRC(mcf,mcDS*mi.pixA,[filtname 'mf.mrc']);
            imwrite(uint8(imscale(rot90(mcf),256,1e-3)),[filtjname 'mf.jpg']);
            imacs(mcf);
        else
            imacs(mc);
        end;
        disp(['Wrote merged images: ' procname 'm.mrc']);
        
        title(mi.baseFilename,'interpreter','none');
        
        % Show the image and spectrum
        set(0,'CurrentFigure',10);  % i.e. figure(10)
        QuickLookSpectrum(mc,mcDS*mi.pixA,[mi.baseFilename 'm.mrc']);
        
        %%
        if removeCrystal
            mcs=meRemoveCrystal(mc,mcDS*mi.pixA,6,1);  % crystal-subtracted image
            WriteMRC(mcs,mcDS*mi.pixA,[procname 'mc.mrc']);
            imwrite(uint8(imscale(rot90(mcs))),[procname 'mc.tif']);
        end;
        %%
        
        if findVesicles
            mi=meFindVesicles(mcs, mi);  % old function, doesn't insert model.
            %             mi=meFindVesicles2(mcs,mi);
            disp('Reconstructing vesicles');
            
            mi=meInsertVesicleModel(mi);
            
            vm=meMakeModelVesicles(mi,mi.imageSize/mcDS);  % downsampled to mc's size
            mcv=mcs-vm;
            
            %%
            subplot(2,3,1);
            imacs(BinImage(mcs,4));
            subplot(2,3,4);
            imacs(BinImage(mcv,4));
            title('Subtracted');
            subplot(2,3,5);
            hist(mi.vesicle.s);
            xlabel('Vesicle amplitude s');
            drawnow;
            WriteMRC(mcv,mcDS*mi.pixA,[procname 'mcv.mrc']);
            imwrite(uint8(imscale(rot90(mcv))),[procname 'mcv.tif']);
            disp(['TIFF pixel size is ' num2str(mi.pixA*mcDS) ' A']);
        end;
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
