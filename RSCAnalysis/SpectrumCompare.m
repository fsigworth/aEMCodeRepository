% SpectrumCompare.m
% Compare 160909 with 191228 images with very similar defocus values
% Or with 200824 images.


% % pathA='/Volumes/Drobo4/siggpu_data/160909/KvLipo121_2/';
pathA='/Volumes/D257/Hideki/160909_spectrAnalysis/';
% % baseName1603='sq02_1_0014_Sep09_19.15.05'; % 3.02 um
% baseNameA='sq02_1_0028_Sep09_19.26.19'; % 2.26um
% infoA=ReadMiFile([pathA 'Info/' baseNameA 'mi.txt']);
% % micA=ReadMRC([path16 'Micrograph_sq02_and_sq81/' baseNameA 'ala.mrc']);
% micA=ReadMRC([pathA 'Micrograph/' baseNameA 'ala.mrc']);
% save AImageData.mat pathA baseNameA infoA micA
% save([pathA 'AImageData.mat'],'pathA', 'baseNameA', 'infoA', 'micA');
% return
% load([pathA 'AImageData.mat']);

%% 
%% 200922 Graphene
        pathB='/Volumes/D257/Hideki/200922/';
        shortPathB='200922';
        cd(pathB);
%         load DefSorted.mat
        load foundMis_broken_graphene
        load AImageData.mat % Data for Image A, but pick this up from pathB!
        
%         jpegDir='jpeg_spect/';
        jpegDir='jpeg_test_spect/'; % test set
        CheckAndMakeDir(jpegDir,1);

suffixes={'m.mrc' 'mv.mrc'};
% for fileIndex=1:20 % graphene set
%     for j=1:2
for fileIndex=1:3 % test set
    for j=1
        % %     pathB='/Volumes/Drobo4/191228.1/';
        % %     shortPathB='191228.1';
        % %     baseNameB='GridSquare_28970644_FoilHole_28972832_02'; % 2.26 um
        
        %%
        %     pathB='/Volumes/D257/Hideki/200816/DefocusNear226/';
        %     shortPathB='200816';
        %     oldPath=pwd;
        %     cd([pathB 'Info/']);
        %     d=dir;
        %     mi=ReadMiFile(d(fileIndex).name);
        %     cd(oldPath);
        %     baseNameB=mi.baseFilename;
        %     % baseNameB='2020-08-17_06_39_41_02505_1_162-5_0002_1-3'; % My favorite image
        
%% using the MotionCor version of file A, to whiten it.
% pathB=pathA;
% shortPathB='motioncorr';
% infoB=infoA;
% micBPath='MotionCorr/job007/movie_frames/';
% baseNameB='Sep09_19_26_19';
% micB=ReadMRC([pathB micBPath baseNameB '.mrc']);
% % micB=micB/mean(micB(:))-1; % normalize
% mic=micB-mean(micB(:));

        %% Get the 200922 graphene comparison variables
        infoB=foundMis{fileIndex};
        baseNameB=infoB.baseFilename;
        mrcNameB=[pathB 'Merged/' baseNameB suffixes{j}];
        disp(mrcNameB);
        
        micB=ReadMRC(mrcNameB);
        
        %%
        mag=infoB.pixA/infoA.pixA;
        % npix=NextNiceNumber(size(mic1902,1)*mag)
        npix=2048;
        mscB=DownsampleGeneral(micB,npix,mag);
        mscA=Crop(micA,npix);
        pixA=infoA.pixA;
        
        figure(1);
        mysubplot(221);
        imags(mscA);
        mysubplot(222);
        imags(mscB);
        
        spA=RadialPowerSpectrum(mscA);
        mysubplot(223);
        semilogy(spA);
        
        spB=RadialPowerSpectrum(mscB);
        mysubplot(224);

        % Define frequencies and fit point xpos, which are where 
%         n=f^2*lambda*def, i.e. zeros of the ctf
        fmax=1/(2*pixA);
        maxN=fmax.^2*22600*.025;
        lambda=EWavelength(200);
        def=infoA.ctf(1).defocus*1e4;
        xpos=sqrt((1:maxN)'/(lambda*def))*npix*pixA;       
        freq=Radius(npix);
        freq1=(0:1023)';
        
        %% Manual fit

%         a1=2;
%         a2=2.5;
%         a3=1.5;
%         a4=3;
%         
%         spBf=11000*spB.*( (1-a1*1e-4*freq1).*(1+freq1.^2*a2*1e-7) ...
%             .*(1+freq1.^3*a3*1e-10).*(1+freq1.^4*a4*1e-13) ).^2;
%         semilogy(min(200,[spBf spA]));
%         hold on;
%         plot(xpos,spBf(round(xpos)),'k.','markersize',10);
%         hold off;
        
        
        %% Simplex fit to minima of spectra
        
        xpts=round(xpos);
        zFreq=freq1(xpts);
        % y19=p(5)*1e4*sp19(xpts);
        y16=spA(xpts);
        scl=[-1e-4 1e-7 1e-10 1e-13];
        p=[1 1 1 1 1];
        p=Simplex('init',p);
        
        for iter=1:1000
            h=1;
            for i=1:4
                h=h.*(1+p(i)*scl(i)*zFreq.^i);
            end;
            yB=h.^2.*spB(xpts)*p(5)*1e4;
            if mod(iter,10)==1
                semilogy(xpos,min(10000,[yB y16]),'-.');
                drawnow;
            end;
            err=(log(yB)-log(y16))./(xpos+400);
            se=err'*err;
            p=Simplex(se);
            title(num2str(p));
        end;
        %% show the final result
        h=1;
        for i=1:4
            h=h.*(1+p(i)*scl(i)*freq.^i);
        end;
        h1=sectr(h);
        yB=(h1.^2.*spB)*p(5)*1e4;
        semilogy(min(10000,[yB spA]),'-');
        drawnow;
        
        %% Having matched the backgrounds, now filter the entire micrograph B to match A
        
        mNewB=real(ifftn(fftn(mscB).*fftshift(h)))*sqrt(p(5)*1e4);
        spBh=RadialPowerSpectrum(mNewB);
        semilogy([spBh spA]);
        
        % The final images are
        % mscA cropped image from A (160909)
        % mNewB resampled image from B
        %%
        fc=.2;
        nSDs=5;
        % mysubplot(221);
        subplot('position',[.03 .3 .46 .68]);
        %     imags(GaussFilt(mscA,fc));
        mAf=GaussFilt(mscA,fc);
        imaga( (mAf-median(mAf(:))) * 128/(nSDs*std(mAf(:))) + 128);
        ylabel(['Pixel size ' num2str(pixA) ' Å']);
        text(50,50,['Pixel size ' num2str(pixA) ' Å'],'color','w','fontsize',12);
        title(['160909/KvLipo121_2 ' baseNameA],'interpreter','none');
        % mysubplot(222);
        subplot('position',[.52 .3 .46 .68]);
        %     imags(GaussFilt(mNewB,fc));
        mBf=GaussFilt(mNewB,fc);
        imaga( (mBf-median(mBf(:))) * 128/(nSDs*std(mBf(:))) + 128);
        
        title([shortPathB '  ' baseNameB suffixes{j}],'interpreter','none');
        % mysubplot(212);
        subplot('position',[.03 .05 .95 .23]);
        semilogy(freq1/(npix*pixA),min(1000,[yB spA]),'-');
        legend(shortPathB,'160909');
        
        jName=[pathB jpegDir baseNameB suffixes{j}(1:2) 'sp.jpg'];
        print('-djpeg',jName);
        
    end;
end;
