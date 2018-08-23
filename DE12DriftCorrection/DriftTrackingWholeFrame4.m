% DriftTrackingWholeFrame4
%  Track drift using the spectral model directly.
% F.S. 4 June 12

% clear;  % Clear out other variables, as this script needs about 4 GB itself.
startIndex=1;
endIndex=inf;
indexStep=1;
postName='';
% % Server
% inDir='/Volumes/cryoEMdata/supervisor/data/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/Tracking/';
% aliDir='/Volumes/cryoEMdata/supervisor/data/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/Aligned/';
% jpegDir='/Volumes/cryoEMdata/supervisor/data/120122/DSR wo carbon 2.5mgml-1 blot1sec/DDD/Aligned/Jpeg/';

% % local disks
% baseDir=AddSlash(uigetdir());
% inDir=[baseDir 'Tracking/'];

inDir='/Volumes/TetraData/EMWork/Anchi/12apr26c/Tracking/';
aliDir=[inDir '../Aligned/'];
jpegDir=[inDir '../Jpeg/'];
postName='wt01';

startIndex=6;
indexStep=10;
endIndex=10;

% 
% % % % For making Figs.1,2 and 4 for the paper use the AMPAR dataset:
% inDir='/Volumes/TetraData/EMWork/Hideki/111229/111229/AMPAR_boxB_slot3/Tracking/';
% aliDir=[inDir '../Aligned/'];
% jpegDir=[inDir '../Jpeg/'];
% startIndex=5  % the set used for the figures
% endIndex=5


% crop=3072;  % For cropping the input images
crop=0;
% crop=2048;

procImages=1;  % Do all the pre-processing
restoreImages=1;  % Align the full-size images
doFullSizeCCs=1;
doComparePowerSpectra=1;  % Show the 'before' and 'after' spectra
saveFigures=1;  % Save copies of figures to the jpeg directory

dsw=4;  % downsampling of working copy
niters=10;
kRelax=1;  % relaxation factor for iterations
cameraIndex=3;  % index for the DE camera
pixA=2.9;  % only needed for indexing crystal spots
showRawCC=0;
useCrystal=0;   % option for streptavidin crystals

%%

names=FindFilenames(inDir,'-s\.mat');
nFiles=numel(names);

if nFiles<1
    error(['Couldn''t find files in ' inDir]);
end;

disp([num2str(nFiles) ' files found.'])



inDir=AddSlash(inDir);
outDir=inDir;

for fileIndex=startIndex:indexStep:min(endIndex,nFiles)
    disp(' ');
    fileIndex
names=FindFilenames(inDir,'-s\.mat');

    inFilename=names{fileIndex};
    disp(inFilename);
    q=regexp(inFilename,'-s');
    if numel(q)<1
        error(['Not a mat file: ' inFilename]);
    end;
    baseFilename=inFilename(1:q-1);  % truncate the end.
    if crop
        baseFilename=[baseFilename 'crop'];
    end;

    
    %% Pre-processing
    
    figure(1);
    SetGrayscale;

    if procImages || restoreImages
        disp('Loading the s struct');
        load([inDir inFilename]);   % Get the images
    end;
    if procImages   % Make working image stacks
        if crop
            s=deCropImageSet(s,crop);  %****************            
        end;
        disp('Fitting the noise model');
        [s f2]=deMakeSensorNoiseModel(s);  % f2 is variables for a figure.
        
        disp('Correcting images');
        [imgs1 imgs2]=deMakeCorrImageSet(s,0,0);  % no prewhitening,
        % no anti-noise
        [nx0 ny0 nim]=size(s.imgs);
        n0=[nx0 ny0];  % Original image size
        disp('Downsampling images');
        nw=n0/dsw;   % working image size
        
        wI1=single(Downsample(imgs1,nw,1));  % downsample the stack
        wI2=single(Downsample(imgs2,nw,1));
        wSpecModel=ifftshift(deEvalSensorNoiseModel(s,dsw));
        %     matches the cross-spectrum ifftn(fftn(i1).*conj(fftn(i2)))/prod(nw)
        
        clear imgs*  % free up memory
        
        wI10=wI1;
        wI20=wI2;
        save([inDir baseFilename '-w'],'wI10','wI20','wSpecModel');  % Save the working image set
    end;
    
    
    
    %%  Initialize for iterations
    
    % load the downsampled image set
    load([inDir baseFilename '-w.mat']);
    
    [nwx nwy nim]=size(wI10);
    nw=[nwx nwy];
    nx0=dsw*nwx;
    ny0=dsw*nwy;
    
    wI1=wI10;  % Working copies of images
    wI2=wI20;
    xsh=zeros(nim,1);
    ysh=zeros(nim,1);
    dxTotal=single(zeros(nim,1));
    dyTotal=single(zeros(nim,1));
    dxAccum=zeros(nim,niters);  % accumulated shifts for plots
    dyAccum=zeros(nim,niters);
    
    fShifts=complex(single(ones([nw nim])));  % Fourier shift arrays
    wISpec=zeros([nw nim]);  % Fourier-shifted spectrum models
    cciLocal=zeros(15,15,niters);  % 1st ccs for figure
    
    %%
    disp('Tracking');
    figure(1);
    for iter=1:niters
        disp(['Iteration ' num2str(iter)]);
        
        % Construct the independent means
        meanInds=1:nim;  % indices of images to sum for means.
        meanN=numel(meanInds)-1;
        wIMean=single(zeros([nw nim]));
        for i=1:nim
            meanIndsX=meanInds(meanInds~=i);  % Delete one entry
            wIMean(:,:,i)=mean(wI2(:,:,meanIndsX),3);
            wISpec(:,:,i)=wSpecModel.*fShifts(:,:,i).*conj(mean(fShifts(:,:,meanIndsX),3));
        end;
        
        %% Compute the weighting filter
        hMask=1;  % Masking is not properly implemented in deGetWeightingFilger
        %     hMask=fuzzymask(nw,2,nw*0.45,nw*.05);
        H1=deGetWeightingFilter2(wI1,wIMean,wISpec,hMask,pixA*dsw,useCrystal);
        %     H1=H1.*fftshift(RadiusNorm(nw));  %%%%
        
        %% show the raw image and the weighting filter.
        subplot(231);
        imacs(mean(wI1,3).*hMask);
        
        subplot(232);
%         xs=(-nwx/4:nwx/4-1)/nwx;
%         ys=(-nwy/4:nwy/4-1)/nwy;
%         contourf(xs,ys,Crop(fftshift(H1),nw/2)');
        xs=(-nwx/2:nwx/2-1)/(nwx*dsw);
        ys=(-nwy/2:nwy/2-1)/(nwy*dsw);
        contourf(xs,ys,fftshift(H1)');
        xlabel('frequency, pix^{-1}');
        title('H1');
        
        %% Compute the shifts from the cross-correlation peaks
        ccSum=zeros(nw);
        for i=1:nim
            crossSpec=H1.*fftn(wI1(:,:,i)).*conj(fftn(wIMean(:,:,i)))/prod(nw);
            cc=fftshift(real(ifftn(crossSpec-H1.*wISpec(:,:,i))));
            ccSum=ccSum+cc;
            [mxv ix iy]=max2di(cc);
            xsh(i)=ix-nwx/2-1;
            ysh(i)=iy-nwy/2-1;
            if i==1
                subplot(236);
                cciLocal(:,:,iter)=Crop(ccSum,size(cciLocal,1));
                imacs(cciLocal(:,:,iter));
                title('Frame 1 CC');
                crossSpec1=fftn(wI1(:,:,i))...
                    .*conj(fftn(wIMean(:,:,i)))/prod(nw)-wISpec(:,:,i);
            end;
            drawnow;
        end;
        % Accumulate the shifts
        dxTotal=dxTotal+kRelax*xsh;  % used for image correction
        dxAccum(:,iter)=dxTotal;     % 2d array for plotting
        dyTotal=dyTotal+kRelax*ysh;
        dyAccum(:,iter)=dyTotal;
        
        %  Display the progress.
        %     Show the averaged cross-correlation, to watch convergence
        if showRawCC  % compute the CC without using the weighting filter.
            for i=1:nim
                crossSpec=fftn(wI1(:,:,i)).*conj(fftn(wIMean(:,:,i)))/prod(nw);
                cc=fftshift(real(ifftn(crossSpec-wISpec(:,:,i))));
                ccSum=ccSum+cc;
            end;
        end;
        %     subplot(236);
        %     imacs(Crop(ccSum,32));
        %     title('Mean cross correlation');
        %     drawnow;
        
        subplot(234);  % show the cumulative translations
        plot(-dsw*dxAccum);
        title('X Shift');
        ylabel('Original pixels');
        xlabel('Frame no.');
        cellTxt=cell(niters,1);
        for j=1:niters
            cellTxt{j}=num2str(j);
        end;
        legend(cellTxt,'location','northeast');
        
        subplot(235);
        plot(-dsw*dyAccum);
        title('Y Shift');
        xlabel('Frame no.');
        
        drawnow;
        
        %% Shift the working images, and update the Fourier shifts
        for i=1:nim
            fSh=FourierShift(nw,-[dxTotal(i) dyTotal(i)]);
            wI1(:,:,i)=real(ifftn(fftn(wI10(:,:,i)).*fSh));
            wI2(:,:,i)=real(ifftn(fftn(wI20(:,:,i)).*fSh));
            fShifts(:,:,i)=fSh;
        end;
    end; % for iter
    %%
    if saveFigures
        set(gcf,'paperpositionmode','auto');
        jName=[jpegDir baseFilename postName '-align.jpg'];
        print('-djpeg','-r300',jName);
    end;
    % put the shift values into the structure.
    shiftX=dsw*dxTotal;
    shiftY=dsw*dyTotal;
    save([inDir baseFilename '-w'],'wI10','wI20','wSpecModel','shiftX','shiftY');  % Save the working image set
    res.shiftX=shiftX;
    res.shiftY=shiftY;
    % print out the shifts
    xShift=dsw*(dxTotal(nim)-dxTotal(1))
    yShift=dsw*(dyTotal(nim)-dyTotal(1))
    
    %%  Restore the full-size images
    
    computeSet2=doFullSizeCCs;
    
    if restoreImages
        figure(1);
        subplot(231);
        if ~exist('s','var')
            disp('Loading the s struct');
            load([inDir inFilename]);   % Get the images
        end;
        disp('Getting full-sized images');
        [startImgs1 startImgs2]=deMakeCorrImageSet(s,0,0);  % not prewhitened
        disp('  shifting them');
        %%
        [nx0 ny0 nim]=size(s.imgs);
        n0=[nx0 ny0];
        imgs10=single(zeros([n0 nim]));
        if doFullSizeCCs
            imgs20=single(zeros([n0 nim]));
        end;
        fshSum=zeros(n0);
        for i=1:nim
            fSh=FourierShift(n0,-dsw*[dxTotal(i) dyTotal(i)]);
            fshSum=fshSum+fSh;
            imgs10(:,:,i)=single(real(ifftn(fftn(double(startImgs1(:,:,i))).*fSh)));
            if doFullSizeCCs  % We'll need the other set of images too.
                imgs20(:,:,i)=single(real(ifftn(fftn(double(startImgs2(:,:,i))).*fSh)));
            end;
            if i==1
                fSh1=fSh;
            end;
        end;
        %%
        % Write the output
        disp(['Writing files ' baseFilename postName '-ali and ' baseFilename postName '-add']);
        WriteMRC(sum(imgs10,3),pixA,[aliDir baseFilename postName '-ali.mrc']);
        disp('Writing imgs10 and imgs20');
        save([inDir baseFilename postName '-ali10.mat'],'imgs10');
        save([inDir baseFilename postName '-ali20.mat'],'imgs20');
        
        WriteJpeg(rot90(sum(imgs10,3)),[jpegDir baseFilename postName '-ali.jpg']);
        WriteJpeg(rot90(sum(startImgs1,3)),[jpegDir baseFilename postName '-add.jpg']);
        disp('done.');
        
        %%  Various diagnostics
        %           First, compare full-size cross-correlations
        if doFullSizeCCs
            disp('Computing full-size cross correlations');
            h0=(fftshift(RadiusNorm(n0)).^.6);
            figure(2);
            SetGrayscale;
            ind=1;
            ndis=[65 49];
            
            meanIndsX=meanInds(meanInds~=ind);  % Delete one entry
            imgs20Mean=mean(imgs20(:,:,meanIndsX),3);
            imgs20Mean=imgs20Mean-mean(imgs20Mean(:));
            startImgs2Mean=mean(startImgs2(:,:,meanIndsX),3);
            startImgs2Mean=startImgs2Mean-mean(startImgs2Mean(:));
            
            sp1=fftn(startImgs1(:,:,ind)).*conj(fftn(startImgs2Mean))/numel(startImgs2Mean);
            sp1(:,1)=sp1(:,1)*.3;
            %         specModel0=ifftshift(deTuneSensorNoiseModel(s,fftshift(sp1)));
            specModel0=ifftshift(deEvalSensorNoiseModel(s,1)*1);
            
            cc1comp=Crop(fftshift(real(ifftn(h0.*(sp1-specModel0)))),ndis);
            subplot(221);  % 1st cc, compensated
            imacs(cc1comp);
            
            cc1ali=Crop(fftshift(real(ifftn(sp1))),ndis);
            subplot(223);  % uncompensated
            imacs(cc1ali);
            
            csp2=fftn(imgs10(:,:,ind)).*conj(fftn(imgs20Mean))/numel(imgs20Mean);
            csp2(:,1)=0;
            sModel0=specModel0.*fSh1.*conj(fshSum-fSh1)/(nim-1);
            ccxc=Crop(fftshift(real(ifftn(h0.*(csp2-sModel0)))),ndis);
            subplot(222);
            imacs(ccxc);
            
            cc1c=Crop(fftshift(real(ifftn(csp2))),ndis);
            subplot(224);
            imacs(cc1c);
            drawnow;
        end;
        
        %%  Compare power spectra
        if doComparePowerSpectra
            subFactor=1.1;  % displayed spectrum looks better if the
            %                             subtracted pixel noise model is scaled up
            %                             slightly.
            
            startImgs1Mean=mean(startImgs1,3);
            startImgs1Mean=startImgs1Mean-mean(startImgs1Mean(:));
            
            endImgs1Mean=mean(imgs10,3);
            endImgs1Mean=endImgs1Mean-mean(endImgs1Mean(:));
            disp('Computing spectra');
            sp1=abs(fftn(double(startImgs1Mean))).^2/prod(n0);
            sp2Raw=abs(fftn(double(endImgs1Mean))).^2/prod(n0);
            modelSpectrum=ifftshift(CCDModelSpectrum2D(3)');
            if any(n0~=[3072 4096])
                modelSpectrum=ifftshift(Downsample(CCDModelSpectrum2D(3)',n0));
            end;
            specModel0=ifftshift(deEvalSensorNoiseModel(s,1));
            fShifts=abs(fshSum/nim).^2;  % shifts effect on spectrum
            sp1=(sp1-specModel0)./modelSpectrum;
            sp2=(sp2Raw-specModel0.*fShifts*subFactor)./modelSpectrum;  % include shift effects.
            
            sp1d=BinImage(fftshift(sp1),4);
            sp2d=BinImage(fftshift(sp2),4);
            sp2RawD=BinImage(fftshift(sp2Raw),4);
            
            res.spectrumAdd=sp1d;
            res.spectrumAli=sp2d;
            %
            %
            pars.spectCropFraction=.8;
            pars.offs1=40;
            pars.offs2=40;
            pars.displayFc=.05;
            pars.satFreq1=.07;
            pars.satFreq2=.07;
            pars.displayExp2=1; %  Spectral value exponent for 2D display
            pars.displayExp1=1;
            pars.lfbin=40;   % number of low-frequency points to ignore in 1d plot
            
            pars.gaussSD=.15;  % Gaussian subtracted from 2d spectrum 1
            pars.gaussAmp=0;
            pars.k1=0;   % exponential change
            pars.fexp1=0.6;  % exponential freq exponent
            
            figure(3);
            SetGrayscale;
            clf;
            disp('Compare power spectra');
            deComparePowerSpectra(sp1d,sp2d,pars);
            if saveFigures
                set(gcf,'paperpositionmode','auto');
                print('-djpeg','-r300',[jpegDir baseFilename postName '-spect.jpg']);
            end;
        end;
        %      deComparePowerSpectra;
        
        % %% Show ACF
        %     ac1=fftshift(real(ifftn(abs(fftn(pimgs)).^2)));
        %     ac2=fftshift(real(ifftn(abs(fftn(sum(imgs10,3))))));
        %     figure(2);
        %     subplot(221);
        %     ccdis=Crop(ac1,32);
        %     ccdis(17,17)=0;
        %     imacs(real(ccdis.^.2))
        %
        
    end; % if restore
    disp('Saving the res structure');
	resFilename=[outDir baseFilename postName '-res']
    save(resFilename,'res');
end; % for