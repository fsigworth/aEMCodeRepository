function mi=k2DriftTrackFcn(mi,pars)

doTracking = pars.forceFitting || (numel(mi.frameShifts)<1);
%%

flagSegments=mi.frameSets;
numDefs=max(1,size(flagSegments,1));
dsw=pars.dsw;
displayOn=1;
nrot=3;  % number of rotations

gainRef=ReadEMFile(mi.gainRefName);

for segIndex=1:numDefs  % defocus stretch
    movName=[mi.moviePath mi.movieFilename];
    disp(['Reading ' movName ' segment ' num2str(segIndex)]);
    %     Get the header information
    [mv, s]=ReadMovie(movName,1,0);
    inds=flagSegments(segIndex,:);
    inds=min(s.nz,inds);
    mi.frameSets(segIndex,:)=inds;  % trim the frame numbers.
    nim=inds(2)-inds(1)+1;
    disp([' frames ' num2str(inds)]);
    pixA=mi.pixA;  % use the pixA value stored in the mi structure
    n=[1 1]*NextNiceNumber(max(s.nx,s.ny));
    
    %     Try to minimize memory use by using singles mostly
    m0=zeros(s.ny,s.nx,nim,'single');
    msk1=true(s.ny,s.nx);
    
    mfMean=0;
    for i=1:nim
        m=single(rot90(ReadMovie(movName,i+inds(1)-1,1),nrot));
        %                 in case there are ADC failures, mask all vertical lines that are zeros
        msk1=repmat(any(m,2),1,s.nx); % blank dead rows in k2 image
        m=m.*msk1;
        pixMean=sum(m(:))/sum(msk1(:));
        %         m1(:,:,i)=Crop((m.*msk0)-pixMean,n);  % pad the entire stack.
        mfMean=mfMean+pixMean;
        m0(:,:,i)=m.*msk1;
    end;
    mfMean=mfMean/nim;
    mi.frameDose=mfMean/(pixA^2);
    mi.imageSize=n;
    m1=zeros([n nim],'single');
    %     We subtract the mean to avoid FFT roundoff errors, and let the
    %     padding be zeros.
    for i=1:nim
        m1(:,:,i)=Crop((m0(:,:,i)-mfMean),n);  % pad the entire stack.
    end;
    
    clear m0  % original size movie
    if doTracking
        
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
        
        % Compute all the ccfs
        fmSum=sum(fm0,3);
        ncc=256;
        
        disp('Computing CCFs');
        sumCs=zeros(n);
        for i=1:nim
            cs=fm0(:,:,i).*conj((fmSum-fm0(:,:,i)));
            sumCs=sumCs+cs;
        end;
        sumCs(1,1:120:n(2))=0;  % Remove the spikes
        ccMean=Crop(fftshift(real(ifftn(sumCs))),ncc)/nim;
        
        % Estimate the artifact amplitude
        
        us=-.21;  % undershoot values
        us2=-.08;
        mul=[1 1 1 1 1; 1 0 0 0 1; 1 1 0 1 1; 1 0 0 0 1; 1 1 1 1 1]';
        pmul=Crop(mul,ncc);
        frac=[us us2 us; 0 1 0; us us2 us]';
        pfrac=Crop(frac,ncc);
        %     Estimate the value of the peak artifact
        ctr=ncc/2+1;
        pkv=ccMean(ctr,ctr)-(pmul(:)'*ccMean(:))/sum(mul(:));
        ccCorr=pkv*pfrac;
        
        % Downsample the images
        disp('Downsampling');
        nw=n/dsw;
        nwx=nw(1);
        nwy=nw(2);
        %     Make the downsampling mask and include the cc artifact compensation
        dsWindow=ifftshift(fuzzymask(nw,2,0.95*nw(1),.1*nw(1)));
        wfCorr=dsWindow.*Cropo(fftn(ifftshift(Crop(ccCorr,n))),nw);
        wF0=zeros([nw nim]);  % Working FTs of images
        wI0=zeros([nw nim]);  % Working means of images
        for i=1:nim
            wF0(:,:,i)=dsWindow.*Cropo(fm0(:,:,i),nw);  % truncated ft
            wI0(:,:,i)=real(ifftn(wF0(:,:,i)));
        end;
        
        clear fm0 % At the end we'll take the FT again, but after gain correction.
        
        %%  Initialize tracking
        wI=wI0;
        
        xsh=zeros(nim,1);
        ysh=zeros(nim,1);
        dxTotal=single(zeros(nim,1));
        dyTotal=single(zeros(nim,1));
        dxAccum=zeros(nim,pars.numIters);  % accumulated shifts for plots
        dyAccum=zeros(nim,pars.numIters);
        
        fShifts=complex(single(ones([nw nim])));  % Fourier shift arrays
        wISpec=zeros([nw nim]);  % Fourier-shifted spectrum models
        %     cciLocal=zeros(15,15,pars.numIters);  % 1st ccs for figure
        
        %
        mdisp(pars.logs,'Tracking');
        for iter=1:pars.numIters
            disp(['Iteration ' num2str(iter)]);
            
            % Construct the independent means
            meanInds=1:nim;  % indices of images to sum for means.
            wIMean=zeros([nw nim]);
            for i=1:nim
                meanIndsX=meanInds(meanInds~=i);  % Delete one entry
                wIMean(:,:,i)=mean(wI(:,:,meanIndsX),3);
                wISpec(:,:,i)=wfCorr.*fShifts(:,:,i).*conj(mean(fShifts(:,:,meanIndsX),3));
            end;
            
            %% Compute the weighting filter
            [H1, rs]=k2GetWeightingFilter3(wI,wIMean,wISpec,5,displayOn);
            
            
            %% Compute the shifts from the cross-correlation peaks
            ccSum=zeros(nw);
            for i=1:nim
                crossSpec=H1.*fftn(wI(:,:,i)).*conj(fftn(wIMean(:,:,i)));
                cc=fftshift(real(ifftn(crossSpec-H1.*wISpec(:,:,i))));
                ccSum=ccSum+cc;
                [mxv ix iy]=max2di(cc);
                xsh(i)=ix-nwx/2-1;
                ysh(i)=iy-nwy/2-1;
            end;
            % Accumulate the shifts
            dxTotal=dxTotal+xsh;  % used for image correction
            dxAccum(:,iter)=dxTotal;     % 2d array for plotting
            dyTotal=dyTotal+ysh;
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
            
            
            trackDisDat=ShowFrameShifts(-dsw*dxAccum, -dsw*dyAccum, inds, iter);
            
            
            if iter>1 && max(abs([xsh; ysh]))<0.1/dsw % shift is less than 0.1 original pixel
                break  % exit from iterations early
            end;
            
            %% Shift the working images, and update the Fourier shifts
            for i=1:nim
                fSh=FourierShift(nw,-[dxTotal(i) dyTotal(i)]);
                wI(:,:,i)=real(ifftn(fftn(wI0(:,:,i)).*fSh));
                fShifts(:,:,i)=fSh;
            end;
            mi.frameShifts{segIndex}=[dxTotal dyTotal]*dsw;
            
        end; % for iter
    else
        trackDisDat=ShowFrameShifts(-mi.frameShifts{segIndex}(:,1),...
            -mi.frameShifts{segIndex}(:,1),inds,1);
        
    end; % doTracking
    %             pause
    %%  Restore the full-size images
    
    if pars.doRestoreImages
        disp('Restoring full-sized images');
        %%
        %         Create the FTs again, but using the gain reference this time.
        if numel(gainRef)>1
            gainRefSq=Crop(gainRef.*msk1,n,0,1);  %%% use msk0 here???
        else
            gainRefSq=Crop(msk1,n);
        end;
        fm0=zeros([n nim],'single');
        %         Restore the mean to m1 for correct gain reference effect, then
        %         subtract it again
        for i=1:nim
            fm0(:,:,i)=fft2(double(((m1(:,:,i)+mfMean).*gainRefSq)-mfMean));
        end;
        
        if pars.doDamageComp
            %             disp('Damage compensation');
            fDamage=k2DamageCompWeights(mi,segIndex);
        else
            fDamage=ones([n nim],'single');
        end;
        fsum=zeros(n);
        fshSum=zeros(n);
        disp('Shifting');
        for i=1:nim
            %             fSh=FourierShift(n,-dsw*[dxTotal(i) dyTotal(i)]);
            fSh=FourierShift(n,-mi.frameShifts{segIndex}(i,:));
            fshSum=fshSum+fSh;
            fsum=fsum+double(fm0(:,:,i)).*fSh.*double(ifftshift(fDamage(:,:,i)));
        end;
        fmSum=sum(fm0,3);
        %%
        ctr=n/2+1;
        imgSum=real(ifftn(fsum))+sum(fDamage(ctr(1),ctr(2),:),3)*mfMean;  % restore the mean
        %         mfMean2=mean(imgSum(:))/nim
        
        % Show the spectra
        spectrumDisplayExp=.1;
        spectrumDisplayLims=[0 5e-3];
        bin=8;
        trackDisDat.img.image=imscale(BinImage(imgSum,4),256,1e-5);
        trackDisDat.img.axis='off';
        trackDisDat.img.title=mi.baseFilename;
        
        pspect0=BinImage(fftshift(abs(fmSum)),bin);
        simg0=Crop(pspect0,n/(bin*2)).^spectrumDisplayExp;
        trackDisDat.sp1.image=imscale(simg0,255,spectrumDisplayLims);
        trackDisDat.sp1.axis='off';
        trackDisDat.sp1.title='Before align: to 1/2 Nyquist';
        
        pspect=BinImage(fftshift(abs(fsum)),bin);
        simg=Crop(pspect,n/(bin*2)).^spectrumDisplayExp;
        trackDisDat.sp2.image=imscale(simg,255,spectrumDisplayLims);
        trackDisDat.sp2.title='After align';
        
        spect1d=Radial(pspect);
        spect1d(1)=NaN;
        trackDisDat.spt.plotx=(0:bin:n(1)/2-1)/(n(1)*pixA);
        trackDisDat.spt.ploty=sqrt(spect1d);
        trackDisDat.spt.yLabel='sqrt(spectrum)';
        trackDisDat.spt.xLabel='frequency, A^{-1}';
        
        
        if displayOn && pars.fig1>0
            figure(1);
            %             DrawFigureFromData(trackDisDat.img,2,3,1);
            DrawFigureFromData(trackDisDat.sp1,2,3,2);
            DrawFigureFromData(trackDisDat.sp2,2,3,3);
            %             DrawFigureFromData(trackDisDat.shiftX,2,3,4);
            %             DrawFigureFromData(trackDisDat.shiftY,2,3,5);
            DrawFigureFromData(trackDisDat.spt,2,3,6);
        end;
        
        %
        %
        %
        %         figure(1);
        %         subplot(2,3,1);
        %         imac(imscale(BinImage(imgSum,4),256,1e-5));
        %         axis off;
        %         title(mi.baseFilename,'interpreter','none');
        %
        %         pspect=BinImage(fftshift(abs(fsum)),bin);
        %         pspect0=BinImage(fftshift(abs(fmSum)),bin);
        %         subplot(2,3,2);
        %         simg=Crop(pspect0,n/(bin*2)).^spectrumDisplayExp;
        %         imac(imscale(simg,255,spectrumDisplayLims));
        %         title('Before align: to 1/2 Nyquist');
        %
        %         subplot(2,3,3);
        %         simg=Crop(pspect,n/(bin*2)).^spectrumDisplayExp;
        %         imac(imscale(simg,255,spectrumDisplayLims));
        %         title('After align');
        %
        %         subplot(2,3,6);
        %         spect1d=Radial(pspect);
        %         spect1d(1)=NaN;
        %         plot((0:bin:n(1)/2-1)/(n(1)*pixA),sqrt(spect1d));
        %         ylabel('sqrt(spectrum)');
        %         xlabel('frequency, A^{-1}');
        %     end;
        %%
        % Write the output
        %             if numDefs>1
        %                 segName=[nameSegment num2str(segIndex)];
        %             else
        %                 segName='';
        %             end;
        outName=[mi.baseFilename pars.nameAli char(pars.nameSegment+segIndex-1)];
        if ~exist(mi.imagePath,'dir')
            mkdir(mi.imagePath);
        end;
        disp(['Writing ' mi.imagePath outName '.mrc']);
        WriteMRC(imgSum,pixA,[mi.imagePath outName '.mrc']);
        WriteJpeg(imgSum,[pars.dirJpeg outName '.jpg']);
        mi.imageFilenames{segIndex}=[outName '.mrc'];
    end;
    if pars.writeGraphics && pars.fig1>0
        set(gcf,'paperpositionmode','auto');
        jName=[pars.dirJpeg outName '-align.jpg'];
        print('-djpeg','-r200',jName);
    end;
    
end;
end

    function ShowFrameShifts(accumShiftsX, accumShiftsY, inds, iter)
        nim=size(accumShiftsX,1);
        ymx=max([-accumShiftsX(:); accumShiftsY(:)]);
        ymn=min([-accumShiftsX(:); accumShiftsY(:)]);
        yspan=(ymx-ymn)*.1;
        
        subplot(234);  % show the cumulative translations
        %                 plot(-dsw*dxAccum(:,1:iter),'.-','markersize',5);
        plot(accumShiftsX(:,1:iter),'.-');
        axis([0 nim+1 ymn-yspan ymx+yspan]);
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
        set(gca,'xtick',1:nim,'xticklabel',xtext,...
            'xgrid','on','ygrid','on');
        cellTxt=cell(iter,1);
        for j=1:iter
            cellTxt{j}=num2str(j);
        end;
        legend(cellTxt,'location','northeast');
        
        subplot(235);
        plot(accumShiftsY(:,1:iter),'.-');
        axis([0 nim+1 ymn-yspan ymx+yspan]);
        title('Y Shift');
        set(gca,'xtick',1:nim,'xticklabel',xtext,...
            'xgrid','on','ygrid','on');
        xlabel('Frame no.');
        
        drawnow;
    end