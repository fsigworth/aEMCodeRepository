%  mode=2;  % For Hideki's standard particles
% mode=3;  % For Hideki's standard particles, 1 image only.
% mode=4;  % For Hideki's sequential files
% mode=5;  % Same, but 1 image only.
mode=7;  % Liguo Leginon data
% mode=10;  % DE12 data
% mode=11;  % Liguo DE12 data
%  mode=12;  % Hideki DE12 AMPAR data
% mode=13;    % Liguo's test data
mode=14;  % Hideki sequential files
% mode=15;  % Hideki DSR aligned images
% mode=15.5 % non-aligned DSR
% % mode=16;  % Qiu-Xing IP3R
% mode=12;
% % mode=15.9;
% mode=14.2  % 120 kV tiff images
MergeAllImages2OldCode
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
            %serverDir='/Volumes/DATA/TEM_DATA/F20/';
            %             imagePath='/Volumes/DATA/TEM_DATA/F20/111006/LND_SoyPC_fr11/';
            serverDir='/Volumes/TetraData/EMWork/Hideki/110628/';
            imagePath=[serverDir 'HS12slot4_GluR2-GFP/'];
            procPath=[serverDir 'Merged/'];
            infoPath=[serverDir 'mi/'];
            serverDir='';
            %procPath=[serverDir imagePath];
            %             procPath=[imagePath];
            %             infoPath=imagePath;
            %             seqMode=0;  % read files of the form 001u4000.dm3 etc.
            %infos=meScanHidekiFiles(serverDir,imagePath,procPath,infoPath,overwrite,seqMode);
            infos=meScanHidekiFiles(serverDir,imagePath,procPath,infoPath,overwrite,0);
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
            initialDefoci=[0.5 4 15];
            %             mainPath='/Volumes/TetraData/EMWork/Liguo/BKfavorite2010/10dec18a/';
            basePath='/Volumes/TetraData/EMWork/Liguo/10oct05a/';
            imagePath=['Images/'];
            procPath=['Merge/'];
            infoPath=['Info/'];
            jpegPath=[procPath 'jpeg/'];
            filtPath=[procPath '../MergeFilt/'];
            filtJPath=[filtPath 'jpeg/'];
            
            infos=meScanLeginonFiles2([basePath imagePath],...
                [basePath procPath], [basePath infoPath],overwrite);
            
        case 10  % DE12 camera data
            numExposures=3;  % double exposure
            initialDefoci=[2 7 15];
            inDir='/Volumes/cryoEMdata/supervisor/data/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/Tracking/';
            
            basePath='/Volumes/TetraData/EMWork/Hideki/';
            imagePath='120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/Images/';
            xPath='120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/';
            procPath=[xPath 'Merge/'];
            infoPath=[xPath 'Info/'];
            infos=meScanDEFiles(numExposures,basePath,imagePath,procPath,infoPath);
            iCamera=3;  % DE12
            ctfOptions.lowB=1;  % lower B factors
            cpe=70;  %%counts per electron??
            overridePixA=2.9
            overrideDose=10
            
        case 11  % LIguo DE12 camera data
            removeCrystal=1;
            numExposures=3;  % double exposure
            initialDefoci=[2 7 15];
            basePath='/Volumes/TetraData/EMWork/Liguo/DE12/ddd-xtal-cryo/Triples/';
            imagePath='Images/';
            procPath=['Merge/'];
            infoPath=['Info/'];
            infos=meScanDEFiles(numExposures,basePath,imagePath,procPath,infoPath);
            iCamera=3;  % DE12
            ctfOptions.lowB=1;  % lower B factors
            cpe=36;  %%counts per electron??
            overridePixA=1.3
            overrideDose=10
            
            
        case 12  % AMPAR DE12 camera data
            numExposures=2;  % double exposure
            initialDefoci=[2 12];
            basePath='/Volumes/TetraData/EMWork/Hideki/111229/111229/AMPAR_boxB_slot3/';
            imagePath='Images/';
            procPath='Merged/';
            infoPath='Info/';
            
            jpegPath=[procPath 'jpeg/'];
            filtPath=[procPath 'filtered/'];
            filtJPath=[filtPath 'jpeg/'];
            
            mi.basePath=basePath;
            mi.imagePath=imagePath;
            mi.procPath=procPath;
            mi.infoPath=infoPath;
            
            %   pick up pairs of successive files:
            infos=meScanDEFiles2(mi,numExposures);
            iCamera=3;  % DE12
            ctfOptions.lowB=1;  % lower B factors
            cpe=70;  %%counts per electron??
            overridePixA=2.9
            overrideDose=10
            startImage=3;
            endIimage=3;
        case 13
            basePath='/Volumes/TetraData/EMWork/Liguo/LW_particle_picker_related_2012_0414/';
            imagePath='Images/';
            procPath='Merged/';
            infoPath='Info/';
            removeCrystal=1;
            initialDefoci=[1 4 15];
            load([basePath '10oct05a_BKcWA18D_00008gr_00020sq_v01_00004hl_v01_00003mi.mat']);
            infos{1}=mi;
            iCamera=1;
            
        case 14
            %             Sequential files, triples
            mi=meCreateMicrographInfoStruct11;
            mi.basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/120428/AMPA-R_lipo_with_carbon/';
            mi.basePath='/Volumes/TetraData/EMWork/Hideki/120428/AMPA-R_lipo_with_carbon/';
            mi.basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/120526/AMPA_R_swelling_washing_after1week/';
            mi.basePath='/Volumes/TetraData/EMWork/Hideki/120711/AMPA_R_dialyzed_centrifuged_sampleA/';

            rexp='\d+_.+\.mrc';  % pick up files such as 004_anytext.mrc
            startImage=2;
            endImage=2;
            singleImage=1;
 
            mi.imagePath='micrograph/';
            mi.infoPath='Info/';
            mi.procPath='Merged/';
            readIfPresent=0;
            infos=meScanHidekiFiles2(mi,3,readIfPresent,rexp);  % 3 exposures
            initialDefoci=[.6 4 15];
            iCamera=1;
            cpe=9;  % SerialEM scaling
            
            basePath=mi.basePath;
            imagePath=mi.imagePath;
            infoPath=mi.infoPath;
            procPath=mi.procPath;
            jpegPath=[mi.procPath 'jpeg/'];
            filtPath=[mi.procPath 'filtered/'];
            filtJPath=[filtPath 'jpeg/'];
           
            
            
            
        case 14.1
            %             Sequential files, triples
            mi=meCreateMicrographInfoStruct11;
            mi.kV=120;
            ctfOptions.kV=120;
            ctfOptions.spFromWholeMicrograph=1;
            mi.basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/120809/hERG/';
            mi.imagePath='micrograph/';
            mi.infoPath='Info/';
            mi.procPath='Merged/';
            readIfPresent=0;
            infos=meScanHidekiFiles2(mi,3,readIfPresent);
            initialDefoci=[1 3 10];
            iCamera=1;
            cpe=9;  % SerialEM scaling
            
            basePath=mi.basePath;
            imagePath=mi.imagePath;
            infoPath=mi.infoPath;
            procPath=mi.procPath;
            jpegPath=[mi.procPath 'jpeg/'];
            filtPath=[mi.procPath 'filtered/'];
            filtJPath=[filtPath 'jpeg/'];
            
            
             case 14.2
            %             Sequential files, triples, tiff files
            mi=meCreateMicrographInfoStruct11;
            mi.basePath='/Volumes/raid3/Hideki/home/SerialEM/F20/SPA/120911/AMPA_R_lipo_PE35x3/';
            mi.imagePath='Micrograph/';
            mi.infoPath='Info/';
            mi.procPath='Merged/';
            readIfPresent=0;
            infos=meScanHidekiFiles2(mi,3,readIfPresent,'.+_\d+\.tif');  % 3 exposures
            initialDefoci=[.6 4 15];
            iCamera=1;
            cpe=9;  % SerialEM scaling
            pixA=2.9;
            mi.kV=120;
            ctfOptions.kV=120;
            ctfOptions.spFromWholeMicrograph=1;
            
            basePath=mi.basePath;
            imagePath=mi.imagePath;
            infoPath=mi.infoPath;
            procPath=mi.procPath;
            jpegPath=[mi.procPath 'jpeg/'];
            filtPath=[mi.procPath 'filtered/'];
            filtJPath=[filtPath 'jpeg/'];
            overridePixA=2.9;
            overrideDose=20;
 
                    case 15.5  % DE12 DSR data without using aligned 1st exposure
            mi=meCreateMicrographInfoStruct11;
            numExposures=3  % single exposure
            initialDefoci=[2 7 15]; %
            mi.basePath='/Volumes/cryoEMdata/supervisor/data/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/';
            
            mi.imagePath='Images/';
            %             aliDir='../Aligned/';
            %             flipFirst=-1  % special case to match drift-tracking
            mi.procPath='MergeNotAligned/';
            mi.infoPath='InfoNotAligned/';
            
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
            
        case 15.9  % DE12 DSR data where we read the existing info files
            mi=meCreateMicrographInfoStruct11;
            mi.basePath='/Volumes/cryoEMdata/supervisor/data/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/';
            mi.infoPath='InfoNotAligned/';
            infos=meLoadInfoFiles([mi.basePath mi.infoPath]);
            mi=infos{1};
            basePath=mi.basePath;
            imagePath=mi.imagePath;
            infoPath=mi.infoPath;
            %             procPath=mi.procPath;
            procPath='MiddleExposure';
            singleImage=2;
            jpegPath=[mi.procPath 'jpeg/'];
            filtPath=[mi.procPath 'filtered/'];
            filtJPath=[filtPath 'jpeg/'];
            iCamera=3;  % DE12
            ctfOptions.lowB=1;  % lower B factors
            ctfOptions.spFromWholeMicrograph=1;
            cpe=70;  %%counts per electron??
            overridePixA=2.9
            overrideDose=10
            %              startImage=11
            mcDS=1  % no downsampling.
            FitAndAlign=0;      % Fit CTF and align.  0 means merge only, assuming that
            overwrite=0;        % That is, should we create new mi files.
            WriteInfo=0;        % Should we store the new mi data
            makeDirectories=0;  % if directories don't exist, create them.
            
            
        case 16  % IP3R data
            mi=meCreateMicrographInfoStruct11;
            numExposures=3  % single exposure
            initialDefoci=[2 4 15]; %
            mi.basePath='/Volumes/TetraData/EMWork/QiuXing/';
            
            mi.imagePath='Images/';
            mi.procPath='Merge/';
            mi.infoPath='Info/';
            
            basePath=mi.basePath;
            imagePath=mi.imagePath;
            infoPath=mi.infoPath;
            procPath=mi.procPath;
            jpegPath=[mi.procPath];
            filtPath=[mi.procPath];
            filtJPath=[filtPath];
            
            rexp='im\d+_.+.tif';   %%% tiff files
            readIfPresent=0;
            infos=meScanHidekiFiles2(mi,3,readIfPresent,rexp);  % 3 exposures
            ctfOptions.lowB=1;  % lower B factors
            ctfOptions.spFromWholeMicrograph=1;
            overridePixA=2.4
            overrideDose=10

