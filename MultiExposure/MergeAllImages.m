% MergeAllImages.m
% Main script for processing raw images, merging them and creating the mi
% files.
% mode=1;  % Liguo's Leginon data
 mode=2;  % For Hideki's standard particles
% mode=3;  % For Hideki's standard particles, 1 image only.
% mode=4;  % For Hideki's sequential files
% mode=5;  % Same, but 1 image only.
% mode=7;  % Liguo's
% mode=100;  % interactive merging

FitAndAlign=1;      % Fit CTF and align.  0 means merge only, assuming that
% alignment parameters are already in mi file.
overwrite=1;        % That is, should we create new mi files.
WriteInfo=1;        % Should we store the new mi data
ignoreOldCTFs=1;    % If a CTF fit already exists, don't use it to set the defocus range.
mcDS=2;             % Composite image downsample ratio
kWiener=.04;        % Wiener filter constant
fc=.1;              % Lowpass filter corner frequency applied to Wiener
useCircMask=0;      % Apply circular antialias filter in meCombineImages
maxNim=inf;         % Limit on number of images to merge
uiPath='';

% load or, if necessary, create the micrograph info structures

startImage=1;       % Use to select some of files in the directory
endImage=inf;

while 1  % repeating for the case of interactive merging.
    
    switch mode
        case 1  % Liguo's RSC data
            removeCrystal=1;
            findVesicles=0;
            initialDefoci=[2 4 15];
            % serverDir='/Volumes/TetraData/EMWork/Liguo/FavoriteBK/';  % main directory where images are stored
            % serverDir='/Volumes/TetraData/EMWork/Liguo/FavoriteBK/04T-BKPo-10dec18a-leginon-part/rawdata/';  % main directory where images are stored
            serverDir='/Volumes/TetraData/EMWork/Liguo/BKfavorite2010/03T-BKPc-10oct02a-Leginon-part/0BKc-good-1002/';
            %         serverDir='';
            imagePath='';         % subdirectory for these images
            infoPath='';            % complete path to where the info files are stored.
            procPath=imagePath;
            infos=meScanLeginonFiles2(serverDir,imagePath,procPath,infoPath,overwrite);
        case 2  % Hideki's standard particles of form <num>u<defocus>.dm3
            removeCrystal=0;
            findVesicles=0;
            ignoreOldCTFs=0;  % We get defocus values from filenames
            initialDefoci=[];  % Set from filename
            serverDir='/Volumes/TetraData/EMWork/Hideki/110628/';
            imagePath='HS12slot4_GluR2-GFP/';
            proucPath=[serverDir imagePath];
            infoPath=imagePath;
            seqMode=0;  % read files of the form 001u4000.dm3 etc.
            infos=meScanHidekiFiles(serverDir,imagePath,procPath,infoPath,overwrite,seqMode);
        case 3 % Hideki's images, but use only the single low-defocus image.
            removeCrystal=0;
            findVesicles=0;
            initialDefoci=[];  % will be obtained from filenames
            serverDir='/Volumes/TetraData/EMWork/Hideki/110628/';
            imagePath='HS12slot4_GluR2-GFP/';
            procPath=[serverDir 'NoMerge/'];
            infoPath='NoMerge/';
            infos=meScanHidekiFiles(serverDir,imagePath,procPath,infoPath,overwrite,0);
            for i=1:numel(infos)  % Change the info to point to only one image file
                infos{i}.rawFile=infos{i}.rawFile(1);
                infos{i}.doses=infos{i}.doses(1);
            end;
        case 4 % Hideki's images that are a simple serial number
            removeCrystal=0;
            findVesicles=0;
            %         initialDefoci=[3 7];  % note that image 25 has a high defocus =7.6,
            initialDefoci=[6 10];  % note that image 25 has a high defocus =7.6,
            %         which is not found by fitctf unless the initialDefocus is set
            %         higher.
            serverDir='/Volumes/TetraData/EMWork/Hideki/110724/';
            imagePath='AR1_1slot3_ARoncarbon/';
            procPath=[serverDir 'Merge/'];
            infoPath='mi/';
            infos=meScanHidekiFiles(serverDir,imagePath,procPath,infoPath,overwrite,1);
        case 5  % Same as 4, but single exposure
            FitAndAlign=1;
            removeCrystal=0;
            findVesicles=0;
            initialDefoci=4;
            serverDir='/Volumes/TetraData/EMWork/Hideki/110724/';
            imagePath='AR1_1slot3_ARoncarbon/';
            procPath=[serverDir 'LowDefocusOnly/'];
            infoPath='LowDefocusOnly/';
            infos=meScanHidekiFiles(serverDir,imagePath,procPath,infoPath,overwrite,1);
            for i=1:numel(infos)  % Change the info to point to only one image file
                infos{i}.rawFile=infos{i}.rawFile(1);
                infos{i}.doses=infos{i}.doses(1);
                if numel(infos{i}.ctf)>1
                    infos{i}.ctf=infos{i}.ctf(1);
                end;
            end;
        case 6 % Hideki's images that are a simple serial number--Novartis LNDs
            removeCrystal=0;
            findVesicles=0;
            initialDefoci=[3 7];  % note that image 25 has a high defocus =7.6,
            %         which is not found by fitctf unless the initialDefocus is set
            %         higher.
            serverDir='/Volumes/TetraData/EMWork/Hideki/110722/';
            imagePath='AR1_1slot1_LND_FULL/';
            procPath=[serverDir 'LNDMerge/'];
            infoPath='LNDmi/';
            infos=meScanHidekiFiles(serverDir,imagePath,procPath,infoPath,overwrite,0);
            
        case 7  % Liguo's RSC data: closed BK
            removeCrystal=1;
            findVesicles=0;
            %         initialDefoci=[1 4 15];
            initialDefoci=[1 4 15];
            initialDefoci=initialDefoci(1:min(numel(initialDefoci),maxNim));
            mainPath='/Volumes/TetraData/EMWork/Liguo/BKfavorite2010/10dec18a/';
            imagePath=[mainPath 'Images/'];
            procPath=[mainPath 'Merge/'];
            infoPath=[mainPath 'Info/'];
            
            
            mainPath='/Volumes/TetraData/EMWork/Liguo/BKfavorite2010/03T-BKPc-10oct02a-Leginon-part/';
            imagePath=[mainPath '0BKc-good-1002/'];
            procPath=[mainPath 'Merge/'];
            infoPath=[mainPath 'Info/'];
            
            infos=meScanLeginonFiles2(imagePath,procPath,infoPath,overwrite);
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
    end;
    
    nmi=numel(infos);
    figure(1);
    figure(10);
    %%
    tic
    for i=startImage:min(nmi,endImage)
        disp(['>>Image index ' num2str(i)]);
        mi=infos{i};  % pick up the structure
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
            fnames{j}=[mi.imagePath mi.imageFilenames{j}];
            %         disp(fnames{j});
            disp(mi.imageFilenames{j});
        end;
        %%
        % Read the images
        [m mi.pixA mi.doses]=meReadImages(fnames);
        sz=size(m);
        mi.nPixels=sz(1:2);
        
        % operate with the pre-whitening filter
        %     disp('mePreWhiten');
        m=mePreWhiten(m);
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
            [rawCTFs mi.spectra]=meFitCTFs(m,mi.pixA,initialDefoci(1:ncts),removeCrystal);
            for j=1:numel(rawCTFs)
                defocus=rawCTFs(j).defocus
            end;
            %%
            % Align the images
            [mi.mergeMatrix mi.ctf]=meAlignMultiExposures(m,mi.pixA,rawCTFs,mi.doses,4);
            %%
            % Store the CTF and alignment information
            if WriteInfo
                save([infoPath mi.baseFilename 'mi.mat'],'mi');
            end;
            
        end;
        %%
        % Create the composite image mc
        [effctf mc mts]=meCombineImages(m,mi.pixA,mi.ctf,mi.mergeMatrix,mi.doses,mcDS,useCircMask);
        mc=-mc;  % go back to reversed contrast
        
        %     [effctf mc]=rsCombineImages(m,mi,mcDS);
        %%
        % Filter the combined image
        mcf=meFilterImage(mc,effctf,mcDS*mi.pixA,kWiener,fc);  % Filter the combined image
        
        % Write out the combined image and filtered version
        procname=[mi.procPath mi.baseFilename];  % Use the processed image directory
        
        WriteMRC(mc,mcDS*mi.pixA,[procname 'm.mrc']);
        imwrite(uint8(imscale(mc)),[procname 'm.jpg']);
        WriteMRC(mcf,mcDS*mi.pixA,[procname 'mf.mrc']);
        imwrite(uint8(imscale(mcf)),[procname 'mf.jpg']);
        imacs(mcf);
        title(mi.baseFilename,'interpreter','none');
        
        % Show the image and spectrum
        set(0,'CurrentFigure',10);  % i.e. figure(10)
        QuickLookSpectrum(mc,mcDS*mi.pixA,[mi.baseFilename 'm.mrc']);
        
        %%
        if removeCrystal
            mcs=meRemoveCrystal(mc,mcDS*mi.pixA,6,1);  % crystal-subtracted image
            WriteMRC(mcs,mcDS*mi.pixA,[procname 'mc.mrc']);
            imwrite(uint8(imscale(mcs)),[procname 'mc.tif']);
        end;
        %%
        
        if findVesicles
mi=meFindVesicles(mcs, mi);
%%
            
            disp('Reconstructing vesicles');
            clf;
            SetGrayscale;
            vm=meMakeModelVesicles(mi,size(mc,1));
            mcv=mcs-vm;
            
            
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
            imwrite(uint8(imscale(mcv)),[procname 'mcv.tif']);
            disp(['TIFF pixel size is ' num2str(mi.pixA*mcDS) ' A']);
        end;
        if WriteInfo
            save([infoPath mi.baseFilename 'mi.mat'],'mi');
        end;
        infos{i}=mi;  % put the info back into the array
        
    end;
    toc
    if mode ~=100
        return
    end;
end;  % while loop

disp([num2str(-startImage+min(nmi,endImage)) ' image sets processed.']);
