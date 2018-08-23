% k2DriftTracker2
% Like k2DriftTracker, but reads info files.
% Create a summed image from a K2 movie, read as an MRC file
% containing raw counts, using a gain reference if available.
% The stored aligned-summed images are padded to 3840 x 3840.
% If there are multiple segments, they are stored separately.
%
%   Various jpegs are stored here in a Jpegs directory 2 levels up.
% experiment/Jpegs/
%
% Both the Micrograph and Jpegs directories are created if they aren't
% already present.
% Derived from DriftTrackingWholeFrame5 for the DE camera.
% F.S. 16 Sep 13
% -modified to mask the black bars in K2 micrographs when ADC channels
% fail.  31 Jan 14.

gNamePost='';
gNameAli='al';
gDamageComp=1;
gNameSegment='a';
gRestoreImages=1;  % Compute the full-size aligned images
gDsw=4;            % Downsampling of working images
gNumIters=4;
gFig1Size=[1000 700];  % Size of main figure, in pixels.

gDirJpeg='Jpeg/';
gSaveFigures=1;

if ~exist('gNames','var') || ~exist('gBatchProcessing','var') || ~gBatchProcessing
    % Have the user select some info files
    [gNames, pathName]=uigetfile('*mi.mat','Select mi files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(gNames)
        gNames={gNames};
    end;
    [rootPath, gDirInfo]=ParsePath(pathName);
    cd(rootPath);
    gStart=1;
else
    gStart=gIndex;
end;

if ~exist(gDirJpeg,'dir')
    mkdir(gDirJpeg);
end;
gNumFiles=numel(gNames);
disp(['Working on ' num2str(gNumFiles) ' files']);
%%
tic
for gIndex=gStart:gNumFiles
    clear -regexp ^[^g] % clear everything that doesn't start with g
    disp([num2str(gIndex) ':  ' gNames{gIndex}]);
    load([gDirInfo gNames{gIndex}]);
    flagSegments=mi.frameSets;
    numDefs=max(1,size(flagSegments,1));
    for segIndex=1:numDefs  % defocus stretch
        %         We seem to have trouble with a memory leak. Hence clear variables
%         clearvars -except f* dir* dsw mi num* gainRef segIndex pa*
        % Set up the main display window
        figure(1);
        SetGrayscale;
        sz=get(0,'screensize');  % Root monitor size.
        gFig1Size=min(gFig1Size,sz(3:4)*0.8);  % Don't allow the figure to be too big.
        pos=sz(3:4)/2-gFig1Size/2;
        set(gcf,'outerposition',[pos gFig1Size]);  % Put the figure in the center of the root monitor.
        [pa, nm]=fileparts(mi.gainRefName);
%                 disp(['Loading the gain reference ' nm]);
        gainRef=ReadEMFile(mi.gainRefName);
        
        movName=[mi.moviePath mi.movieFilename];
        disp(['Reading ' movName ' segment ' num2str(segIndex)]);
        %     Get the header information
        [mv, s]=ReadMovie(movName,1,0);
        inds=flagSegments(segIndex,:);
        inds=min(s.nz,inds);
        nim=inds(2)-inds(1)+1;
        disp([' frames ' num2str(inds)]);
        pixA=mi.pixA;  % use the pixA value stored in the mi structure
        
        %     Try to minimize memory use by using singles mostly
        m0=zeros(s.ny,s.nx,nim,'single');
        msk1=true(s.ny,s.nx);
        
        mfMean=0;
        for i=1:nim
            m=single(rot90(ReadMovie(movName,i+inds(1)-1,1),3));
            %                 in case there are ADC failures, mask all vertical lines that are zeros
            msk0=repmat(any(m,2),1,s.nx);
            msk1=msk1&msk0;
            m=m.*msk0;
            pixMean=sum(m(:))/sum(msk0(:));
            %                 m(~msk0)=pixMean;
            mfMean=mfMean+pixMean;
            m0(:,:,i)=m.*msk0;
        end;
        mfMean=mfMean/nim;
        %     dn=(n(2)-n0(2))/2;  % First element shift
        n=[1 1]*NextNiceNumber(max(s.nx,s.ny));
        mi.imageSize=n;
        m1=zeros([n nim],'single');
        for i=1:nim
            m1(:,:,i)=Crop((m0(:,:,i)-mfMean),n);  % pad the entire stack.
        end;
        clear m0
        
        % Compute all the fts
        disp('Computing FTs');
        fm0=zeros([n nim],'single');
        for i=1:nim
            fm0(:,:,i)=fft2(double(m1(:,:,i)));
        end;
        %%
        clf;
        subplot(2,3,1);
        imac(imscale(BinImage(sum(m1,3),8),256,1e-5));
        axis off;
        title(mi.baseFilename,'interpreter','none');
        clear m1
        
        % Compute all the ccfs
        fmSum=sum(fm0,3);
        crossSpect=zeros([n nim],'single');
        ncc=256;
        ccs0=zeros(ncc,ncc,nim,'single');
        ccs =zeros(ncc,ncc,nim,'single');
        ctr=ncc/2+1;
        disp('Computing CCFs');
        subplot(2,3,2);
        for i=1:nim
            cs=fm0(:,:,i).*conj((fmSum-fm0(:,:,i)));
            cs(1,1:120:n(2))=0;  % Remove the spikes
            crossSpect(:,:,i)=cs;
            cc=Crop(fftshift(real(ifftn(double(cs)))),ncc);
            ccs0(:,:,i)=cc;
            imacs(GaussFilt(Crop(cc,64),.2));
            title(i);
            drawnow;
        end;
        
        % Estimate the artifact amplitude
        
        us=-.21;  % undershoot values
        us2=-.08;
        mul=[1 1 1 1 1; 1 0 0 0 1; 1 1 0 1 1; 1 0 0 0 1; 1 1 1 1 1]';
        pmul=Crop(mul,ncc);
        frac=[us us2 us; 0 1 0; us us2 us]';
        pfrac=Crop(frac,ncc);
        amps=zeros(nim,1);
        ccMean=mean(ccs0,3);
        %     Estimate the value of the peak artifact
        pkv=ccMean(ctr,ctr)-(pmul(:)'*ccMean(:))/sum(mul(:));
        ccCorr=pkv*pfrac;
        subplot(2,3,3);
        for i=1:nim
            ccs(:,:,i)=ccs0(:,:,i)-ccCorr;
            imacs(GaussFilt(Crop(ccs(:,:,i),64),.2));
            title('Corrected');
            drawnow;
        end;
        %%
        
        % Downsample the images
        disp('Downsampling');
        nw=n/gDsw;
        nwx=nw(1);
        nwy=nw(2);
        dsWindow=ifftshift(fuzzymask(nw,2,0.95*nw(1),.1*nw(1)));
        wfCorr=dsWindow.*Cropo(fftn(ifftshift(Crop(ccCorr,n))),nw);
        wF0=zeros([nw nim]);
        wI0=zeros([nw nim]);
        for i=1:nim
            wF0(:,:,i)=dsWindow.*Cropo(fm0(:,:,i),nw);  % truncated ft
            wI0(:,:,i)=real(ifftn(wF0(:,:,i)));
        end;
        
        %%  Initialize tracking
        wI=wI0;
        
        xsh=zeros(nim,1);
        ysh=zeros(nim,1);
        dxTotal=single(zeros(nim,1));
        dyTotal=single(zeros(nim,1));
        dxAccum=zeros(nim,gNumIters);  % accumulated shifts for plots
        dyAccum=zeros(nim,gNumIters);
        
        fShifts=complex(single(ones([nw nim])));  % Fourier shift arrays
        wISpec=zeros([nw nim]);  % Fourier-shifted spectrum models
        cciLocal=zeros(15,15,gNumIters);  % 1st ccs for figure
        
        %
        disp('Tracking');
        for iter=1:gNumIters
            disp(['Iteration ' num2str(iter)]);
            
            % Construct the independent means
            meanInds=1:nim;  % indices of images to sum for means.
            meanN=numel(meanInds)-1;
            wIMean=zeros([nw nim]);
            for i=1:nim
                meanIndsX=meanInds(meanInds~=i);  % Delete one entry
                wIMean(:,:,i)=mean(wI(:,:,meanIndsX),3);
                wISpec(:,:,i)=wfCorr.*fShifts(:,:,i).*conj(mean(fShifts(:,:,meanIndsX),3));
            end;
            
            %% Compute the weighting filter
            [H1, rs]=k2GetWeightingFilter3(wI,wIMean,wISpec,5);
            
            %% show the weighting filter.
            %         subplot(233);
            %         plot(Radial2(fftshift(H1)));
            %         title('H1');
            
            %% Compute the shifts from the cross-correlation peaks
            ccSum=zeros(nw);
            for i=1:nim
                crossSpec=H1.*fftn(wI(:,:,i)).*conj(fftn(wIMean(:,:,i)));
                cc=fftshift(real(ifftn(crossSpec-H1.*wISpec(:,:,i))));
                ccSum=ccSum+cc;
                [mxv ix iy]=max2di(cc);
                xsh(i)=ix-nwx/2-1;
                ysh(i)=iy-nwy/2-1;
                %             if i==1
                %                 subplot(236);
                %                 cciLocal(:,:,iter)=Crop(ccSum,size(cciLocal,1));
                %                 imacs(cciLocal(:,:,iter));
                %                 title('Frame 1 CC');
                %                 crossSpec1=fftn(wI(:,:,i))...
                %                     .*conj(fftn(wIMean(:,:,i)))-wISpec(:,:,i);
                %             end;
                %             drawnow;
            end;
            kRelax=1;
            % Accumulate the shifts
            dxTotal=dxTotal+kRelax*xsh;  % used for image correction
            dxAccum(:,iter)=dxTotal;     % 2d array for plotting
            dyTotal=dyTotal+kRelax*ysh;
            dyAccum(:,iter)=dyTotal;
            %  Display the progress.
            %     Show the averaged cross-correlation, to watch convergence
            showRawCC=0;
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
            
            ymx=max(-gDsw*[dxAccum(:);dyAccum(:)]);
            ymn=min(-gDsw*[dxAccum(:);dyAccum(:)]);
            yspan=(ymx-ymn)*.1;
            
            subplot(234);  % show the cumulative translations
            %                 plot(-dsw*dxAccum(:,1:iter),'.-','markersize',5);
            plot(-gDsw*dxAccum(:,1:iter),'.-');
            axis([0 nim+1 ymn-yspan ymx+yspan]);
            xvals=inds(1):inds(2);
            xtext=cell(1);
            if nim<5
                modVal=1;
            else
                modVal=5;
            end;
            for k=inds(1):inds(2)
                k1=k-inds(1)+1;
                if mod(k,modVal)==0
                    xtext{k1}=num2str(k);
                else
                    xtext{k1}='';
                end;
            end;
            title('X Shift');
            ylabel('Original pixels');
            xlabel('Frame no.');
            set(gca,'xtick',1:size(dxAccum,1),'xticklabel',xtext,...
                'xgrid','on','ygrid','on');
            cellTxt=cell(iter,1);
            for j=1:iter
                cellTxt{j}=num2str(j);
            end;
            legend(cellTxt,'location','northeast');
            
            subplot(235);
            plot(-gDsw*dyAccum(:,1:iter),'.-');
            axis([0 nim+1 ymn-yspan ymx+yspan]);
            title('Y Shift');
            set(gca,'xtick',1:size(dxAccum,1),'xticklabel',xtext,...
                'xgrid','on','ygrid','on');
            xlabel('Frame no.');
            
            drawnow;
            
            if iter>1 && max(abs([xsh; ysh]))<0.1/gDsw % shift is less than 0.1 original pixel
                break  % exit from iterations early
            end;
            
            %% Shift the working images, and update the Fourier shifts
            for i=1:nim
                fSh=FourierShift(nw,-[dxTotal(i) dyTotal(i)]);
                wI(:,:,i)=real(ifftn(fftn(wI0(:,:,i)).*fSh));
                fShifts(:,:,i)=fSh;
            end;
        end; % for iter
        %             pause
        %%  Restore the full-size images
        
        if gRestoreImages
            disp('Getting full-sized images');
            %%
            if gDamageComp
                disp('Damage compensation');
                fDamage=k2DamageCompWeights(mi,segIndex);
            else
                fDamage=ones([n nim],'single');
            end;
            fsum=zeros(n);
            fshSum=zeros(n);
            disp('Shifting');
            for i=1:nim
                fSh=FourierShift(n,-gDsw*[dxTotal(i) dyTotal(i)]);
                fshSum=fshSum+fSh;
                fsum=fsum+double(fm0(:,:,i)).*fSh.*double(ifftshift(fDamage(:,:,i)));
            end;
            %%
            if numel(gainRef)>1
                gainRefSq=Crop(gainRef.*msk1,n,0,1);  %%% use msk0 here???
            else
                gainRefSq=msk1;
            end;
            ctr=n/2+1;
            imgSum=gainRefSq.*real(ifftn(fsum))+sum(fDamage(ctr(1),ctr(2),:),3)*mfMean;  % restore the mean

            % Show the spectra
            spectrumDisplayExp=.1;
            spectrumDisplayLims=[0 5e-3];
            bin=8;
            figure(1);
            subplot(2,3,1);
            imac(imscale(BinImage(imgSum,4),256,1e-5));
            axis off;
            title(mi.baseFilename,'interpreter','none');
            
            pspect=BinImage(fftshift(abs(fsum)),bin);
            pspect0=BinImage(fftshift(abs(fmSum)),bin);
            subplot(2,3,2);
            simg=Crop(pspect0,n/(bin*2)).^spectrumDisplayExp;
            imac(imscale(simg,255,spectrumDisplayLims));
            title('Before align: to 1/2 Nyquist');
            
            subplot(2,3,3);
            simg=Crop(pspect,n/(bin*2)).^spectrumDisplayExp;
            imac(imscale(simg,255,spectrumDisplayLims));
            title('After align');
            
            subplot(2,3,6);
            spect1d=Radial(pspect);
            spect1d(1)=NaN;
            plot((0:bin:n(1)/2-1)/(n(1)*pixA),sqrt(spect1d));
            ylabel('sqrt(spectrum)');
            xlabel('frequency, A^{-1}');
            
            %%
            % Write the output
%             if numDefs>1
%                 segName=[nameSegment num2str(segIndex)];
%             else
%                 segName='';
%             end;
            outName=[mi.baseFilename gNameAli char(gNameSegment+segIndex-1)];
            if ~exist(mi.imagePath)
                mkdir(mi.imagePath);
            end;
            disp(['Writing ' mi.imagePath outName '.mrc']);
            WriteMRC(imgSum,pixA,[mi.imagePath outName '.mrc']);
            WriteJpeg(rot90(imgSum),[gDirJpeg outName '.jpg']);
            mi.imageFilenames{segIndex}=[outName '.mrc'];
            if gSaveFigures
                set(gcf,'paperpositionmode','auto');
                jName=[gDirJpeg outName '-align.jpg'];
                print('-djpeg','-r200',jName);
            end;
            shiftX=dxTotal*gDsw;
            shiftY=dyTotal*gDsw;
            mi.frameShifts{segIndex}=[dxTotal dyTotal]*gDsw;
            mi.pixA=pixA;
            mi.frameDose=pixMean/(pixA^2);
            disp(['Updating ' gNames{gIndex}]);
            save([gDirInfo gNames{gIndex}],'mi');
            disp('----');
            % % %             save([AddSlash(dirMovie) baseFilename '-shifts' segName '.mat'],'shiftX','shiftY');  % Save the working image set
        end;
        %%
        %
    end;
end;
toc;
gIndex=1;