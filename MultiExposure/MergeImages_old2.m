% MergeImages2.m
% Same as the old MergeImages, but operates on already extant mi files

WriteInfo=1;        % Store the new mi file
makeDirectories=1;  % if directories don't exist, create them.
ignoreOldCTFs=1;    % If a CTF fit already exists, don't use it to set the defocus range.
fitAndAlign=1;      % Fit CTF and align images before merging.
mcDS=2;             % Composite image downsample ratio
useCircMask=1;      % Apply circular antialias filter in meCombineImages
maxNim=inf;         % number of image sets to merge
initialDefoci=[1 10];

startImage=1;
endImage=inf;

ctfOptions=struct;
ctfOptions.lowB=1;  % lower B factors
ctfOptions.spFromWholeMicrograph=1;

if ~exist('gNames','var') || ~exist('gBatchProcessing','var') || ~gBatchProcessing
        % Select info files from the file selector

        [gNames, pathName]=uigetfile('*mi.mat','Select info files','multiselect','on');
        if isnumeric(pathName) % File selection was cancelled
            return
        end;
        if ~iscell(gNames)
            gNames={gNames};
        end;
        [rootPath, gInfoPath]=ParsePath(pathName);
        cd(rootPath);
else
    rootPath=pwd;
end;                

%% Main program %%

figure(10);
SetGrayscale;
figure(1);
SetGrayscale;
%%

nmi=numel(gNames);
lastImage=min(nmi,endImage);

tic
disp(['Merging image sets ' num2str(startImage) ' to ' num2str(lastImage)]);
for imIndex=startImage:lastImage
    
    disp(['>>Image index ' num2str(imIndex)]);
    load([gInfoPath gNames{imIndex}]);  
    mi.basePath=rootPath;
    
    if imIndex==1  % set up the directories
        jpegPath=[mi.procPath 'jpeg/'];
        filtPath=[mi.procPath 'filtered/'];
        filtJPath=[filtPath 'jpeg/'];
        
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
        
        modelSpectrum=CCDModelSpectrum2D(mi.camera);  % handle DE-12 or CCD
        [nCx, nCy]=size(modelSpectrum);  % Use the non-normalized spectrum, as
        % mePreWhiten normalizes it.
        
    end;
    
    
    
    %%
    % Read the images
    [m, mi.pixA, mi.doses]=meReadImagesNorm(mi,mi.cpe);
    nim=size(m,3);
    mi.weights=single(ones(1,nim));
    doses=mi.doses
    pixA=mi.pixA
    ctfOptions.kV=mi.kV;
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
        
%         if doRemoveGoldParticles
%             m2=RemoveGoldParticles(m);
%         else
            m2=m;
%         end;
        
        for j=nim:-1:1;
            [rawCTF]=meFitCTF(m2(:,:,j),mi.pixA,...
                initialDefoci(j),0,ctfOptions);  % Actual CTF fitting
            if j==nim
                rawCTFs=rawCTF;
            end;
                rawCTFs(j)=rawCTF;
            defocus=[rawCTFs(j).defocus rawCTFs(j).deltadef]
            jctfName=[mi.basePath filtJPath mi.baseFilename num2str(j) '-spect.jpg'];
            print('-djpeg','-r150',jctfName);  % save the CTF window.
        end;
        
        %%
        % Align the images
        % use downsampling by 8!
        [mi.mergeMatrix, mi.ctf, Ps]=meAlignMultiExposures(m,mi.pixA,rawCTFs,mi.doses,16);
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
    [effctf, mc, mts]=meCombineImages2(m,mi,mcDS,useCircMask);
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
%     if removeCrystal
%         mcs=meRemoveCrystal(mc,mcDS*mi.pixA,6,1);  % crystal-subtracted image
%         WriteMRC(mcs,mcDS*mi.pixA,[procname 'mc.mrc']);
%         WriteJpeg(mcs,[procname 'mc.jpg']);
%     end;
    
    %%
    if WriteInfo
        infoFullName=meSaveMiFile(mi);
        disp(['Info file written: ' infoFullName]);
    end;
    disp(' ');
    
end;
toc

disp([num2str(lastImage-startImage+1) ' image sets processed.']);
disp(' ');
