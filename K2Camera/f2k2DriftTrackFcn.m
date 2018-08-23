function mi=f2k2DriftTrackFcn(mi,pars)


doTracking = pars.forceTracking || ~isfield(mi,'frameShifts') || numel(mi.frameShifts)<1;

flagSegments=mi.frameSets;
numSegments=max(1,size(flagSegments,1));
dsw=pars.dsw;
displayOn=pars.showGraphics;
showRawCC=1;
disDat=struct;
if pars.ggDamageModel
    if mi.kV>200
        mi.damageModelCode='.245*f.^-1.665+2.81; % Grant&Grigorieff 300 kV';
    else
        mi.damageModelCode='.184*f.^-1.665+2.1; % Grant&Grigorieff 200 kV';
    end;
end;

if pars.k2Mode
    if pars.oldK2==1
        nrot=3;
    elseif pars.oldK2==0
        nrot=2;
    else
        nrot=0;
    end;
else
    nrot=0;  % number of 90-degree rotations
end;

ok=true;
for segIndex=1:numSegments  % index of defocus stretch
    if isa(mi.movieFilename,'cell')
        movName=[mi.moviePath mi.movieFilename{segIndex}];
    else
        movName=[mi.moviePath mi.movieFilename];
    end;
    mdisp(pars.logs,['Reading ' movName ' segment ' num2str(segIndex)]);
    %     Get the header information
    [name,ok]=CheckForImageOrZTiff(movName);
    if ~ok
        disp(['--file not found.  pwd= ' pwd]);
        return
    end;
    disp(['Reading ' name]);
    [~, s]=ReadMovie(name,1,0);
    %     Figure out which frames to read
    inds=flagSegments(segIndex,:);
    inds=min(s.nz,inds);
    inds1=inds; % frames to search for dose calculations
    if segIndex<numSegments
        inds1(2)=max(inds(2),flagSegments(segIndex+1,1)-1); % look ahead in single movie
    end;
    inds1=min(s.nz,inds1);
    
    mi.frameSets(segIndex,:)=inds;  % trim the frame numbers.
    nim=inds(2)-inds(1)+1;     % no. of frames to process
    nim1=inds1(2)-inds1(1)+1;  % no. of frames for computing doses
    mdisp(pars.logs,[' frames ' num2str(inds1)]);
    pixA=mi.pixA;  % use the pixA value stored in the mi structure
    n=[1 1]*NextNiceNumber(max(s.nx,s.ny));
    
    m1=zeros([n nim1],'single');
    frameMeans=zeros(nim1,1);
    
    %     gMsk=rot90(false(s.nx,s.ny),nrot);  % this will be the global mask
    if pars.k2Mode % subtract the mean and pad to size n
        for i=1:nim1
            [m,~,ok]=ReadMovie(name,i+inds1(1)-1,1);
            if ~ok
                return
            end;
            m=single(rot90(m,nrot));
            msk1=repmat(any(m,2),1,size(m,2)); % blank dead rows in k2 image
            %             (note that blanking makes an assumption about rotations)
            m=m.*msk1;
            frameMeans(i)=sum(m(:))/sum(msk1(:));
            if i==1
                sigma=std(m(:));
            end;
            lMsk=m>(frameMeans(i)+8*sigma);  % local mask
            %             gMsk=gMsk | lMsk;  % accumulate masked points.
            m(lMsk)=frameMeans(i);  % remove outliers
            m1(:,:,i)=Crop(m-frameMeans(i),n);
        end;
    else
        for i=1:nim1  % simply subtract the mean
            [m,~,ok]=ReadMovie(name,i+inds1(1)-1,1);
            frameMeans(i)=mean(m(:));
        end;
        m1(:,:,i)=(m-frameMeans(i));  % normalize the frame
    end;
    
    mi.frameDose(inds(1):inds(1)+nim1-1)=frameMeans/(mi.pixA^2)/mi.cpe;
    mi.imageSize=n;
    
    if doTracking
        
        %         mdisp(pars.logs,'Computing FTs');
        fm0=zeros([n nim],'single');
        for i=1:nim
            fm0(:,:,i)=fft2(double(m1(:,:,i)/frameMeans(i)));  % normalize by frame dose
        end;
        if pars.k2Mode
            fm0(:,1,:)=0;  % zero the f_y=0 component to remove horizontal artifact
        end;
        %       Graphics: first entry into trackDisDat
        disDat.img.image=imscale(BinImage(sum(m1,3),8),256,1e-5);
        disDat.img.axis='off';
        disDat.img.titleTex=mi.baseFilename;
        
        if displayOn
            clf;
            subplot(2,3,1);
            imaga(disDat.img.image);
            axis off;
            title(mi.baseFilename,'interpreter','none');
        end;
        
        % Compute all the ccfs
        fmSum=sum(fm0,3);
        ncc=256;
        
        %         mdisp(pars.logs,'Computing CCFs');
        sumCs=zeros(n);
        for i=1:nim
            cs=fm0(:,:,i).*conj((fmSum-fm0(:,:,i)));
            sumCs=sumCs+cs;
        end;
        
        %         Prepare for downsampling
        nw=n/dsw;
        nwx=nw(1);
        nwy=nw(2);
        dsWindow=ifftshift(fuzzymask(nw,2,0.95*nw(1),.1*nw(1)));
        
        % CC artifact correction
        if pars.k2Mode
            ctr=ncc/2+1;

            switch pars.ccCorrection
                case 'k2'
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
            pkv=ccMean(ctr,ctr)-(pmul(:)'*ccMean(:))/sum(mul(:));
            ccCorr=pkv*pfrac;
            wfCorr=dsWindow.*Cropo(fftn(ifftshift(Crop(ccCorr,n))),nw);
                case 'delta'
                    mul=Crop(ones(3,3)/8,ncc);
                                ccMean=Crop(fftshift(real(ifftn(sumCs))),ncc)/nim;
                    pkv=ccMean(ctr,ctr);
                    ccCorr=pkv-mul(:)'*ccMean(:);
                    wfCorr=dsWindow*ccCorr;
            end;
        else
            wfCorr=0;  % Do no correction
        end;
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
            %             mdisp(pars.logs,['Iteration ' num2str(iter)]);
            
            % Construct the independent means
            meanInds=1:nim;  % indices of images to sum for means.
            wIMean=zeros([nw nim]);
            for i=1:nim
                meanIndsX=meanInds(meanInds~=i);  % Delete one entry
                wIMean(:,:,i)=mean(wI(:,:,meanIndsX),3);
                if pars.k2Mode  % compute artifact correction
                    wISpec(:,:,i)=wfCorr.*fShifts(:,:,i).*conj(mean(fShifts(:,:,meanIndsX),3));
                end;
            end;
            
            %% Compute the weighting filter
            [H1, rs]=k2GetWeightingFilter3(wI,wIMean,wISpec,5,displayOn);
            
            
            %% Compute the shifts from the cross-correlation peaks
            ccSum=zeros(nw);
            for i=1:nim
                crossSpec=H1.*fftn(wI(:,:,i)).*conj(fftn(wIMean(:,:,i)));
                cc=fftshift(real(ifftn(crossSpec-H1.*wISpec(:,:,i))));
                ccSum=ccSum+cc;
  imags(Crop(cc,32)); title(i); pause(0.2);
                [~, ix, iy]=max2di(cc);
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
            if showRawCC  % compute the CC without using the weighting filter.
                for i=1:nim
                    crossSpec=fftn(wI(:,:,i)).*conj(fftn(wIMean(:,:,i)))/prod(nw);
                    cc=fftshift(real(ifftn(crossSpec-wISpec(:,:,i))));
                    ccSum=ccSum+cc;
                end;
            end;
            
            
            disDat=ShowFrameShifts(-dsw*dxAccum, -dsw*dyAccum, inds, iter, displayOn,disDat);
            
            
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
        disp('Using the existing shift information.');
        disDat=ShowFrameShifts(-mi.frameShifts{segIndex}(:,1),...
            -mi.frameShifts{segIndex}(:,1),inds,1,displayOn,disDat);
        
    end; % doTracking
    %             pause
    %%  Restore the full-size images
    
    if pars.doRestoreImages
        mdisp(pars.logs,'Restoring full-sized images');
        %%
        %         Create the FTs again, but using the gain reference this time.
        fm0=zeros([n nim],'single');
        if pars.k2Mode
            if exist(mi.gainRefName,'file')
                gRef=ReadEMFile(mi.gainRefName);
                %                 gRef(gMsk)=1;  % use the global mask
                gRef=Crop(gRef,n,0,1);
                gainRef=RemoveOutliers(Crop(gRef,n,0,1));
            else
                gainRef=1;
            end;
            gainRefSq=gainRef.*Crop(msk1,n,0,1);  %%% use msk0 here???
            for i=1:nim
                fm0(:,:,i)=fft2(double(((m1(:,:,i)+frameMeans(i))...
                    .*gainRefSq)-frameMeans(i)));
            end;
            
        else  % F2 images are already gain-corrected
            for i=1:nim
                fm0(:,:,i)=fft2(double(m1(:,:,i)));
            end;
        end;
        
        if pars.doDamageComp
            %             mdisp(pars.logs,'Damage compensation');
            fDamage=k2DamageCompWeights(mi,segIndex);
        else
            fDamage=ones([n nim],'single');
        end;
        fsum=zeros(n);
        fshSum=zeros(n);
        if pars.writeStack
           rFrames=ones([n nim],'single');
        end;
        %         mdisp(pars.logs,'Shifting');
        for i=1:nim
            %             fSh=FourierShift(n,-dsw*[dxTotal(i) dyTotal(i)]);
            fSh=FourierShift(n,-mi.frameShifts{segIndex}(i,:));
            fshSum=fshSum+fSh;
            fsum=fsum+double(fm0(:,:,i)).*fSh.*double(ifftshift(fDamage(:,:,i)));
            if pars.writeStack % compute the shifted frames (no damage comp)
                rFrames(:,:,i)=rot90(real(ifftn(fm0(:,:,i).*fSh))+frameMeans(i),pars.finalRot90);
            end;
        end;
        fmSum=sum(fm0,3);
        %%  Restore the mean value and perform a final rotation
        ctr=n/2+1;
        dcDamages=fDamage(ctr(1),ctr(2),:);
        imgSum=real(ifftn(fsum))+frameMeans(1:nim)'*dcDamages(:);  % restore the mean
        if pars.finalRot90>0
            imgSum=rot90(imgSum,pars.finalRot90);
            mi.finalRot90=pars.finalRot90;
        end;
        
        %% ----display code: show the spectra------
        spectrumDisplayExp=.1;
        spectrumDisplayLims=[0 5e-3];
        bin=8;
        disDat.img.image=imscale(BinImage(imgSum,4),256,1e-5);
        disDat.img.axis='off';
        disDat.img.title=mi.baseFilename;
        
        pspect0=BinImage(fftshift(abs(fmSum)),bin);
        simg0=Crop(pspect0,n/(bin*2)).^spectrumDisplayExp;
        disDat.sp1.image=imscale(simg0,255,spectrumDisplayLims);
        disDat.sp1.axis='off';
        disDat.sp1.title='Before align: to 1/2 Nyquist';
        
        pspect=BinImage(fftshift(abs(fsum)),bin);
        simg=Crop(pspect,n/(bin*2)).^spectrumDisplayExp;
        disDat.sp2.image=imscale(simg,255,spectrumDisplayLims);
        disDat.sp2.title='After align';
        
        spect1d0=Radial(pspect0);
        spect1d0(1)=NaN;
        spect1d=Radial(pspect);
        spect1d(1)=NaN;
        disDat.spt.plotx=(0:bin:n(1)/2-1)/(n(1)*pixA);
        disDat.spt.ploty=sqrt([spect1d0 spect1d]);
        disDat.spt.yLabel='sqrt(spectrum)';
        disDat.spt.xLabel='frequency, A^{-1} (red: aligned)';
        
        if displayOn && pars.fig1>0
            figure(1);
            %             DrawFigureFromData(trackDisDat.img,2,3,1);
            DrawFigureFromData(disDat.sp1,2,3,2);
            DrawFigureFromData(disDat.sp2,2,3,3);
            %             DrawFigureFromData(trackDisDat.shiftX,2,3,4);
            %             DrawFigureFromData(trackDisDat.shiftY,2,3,5);
            DrawFigureFromData(disDat.spt,2,3,6);
            %
            %             subplot(2,3,1);
            %             imac(imscale(BinImage(imgSum,4),256,1e-5));
            %             axis off;
            %             title(mi.baseFilename,'interpreter','none');
            %
            % %             pspect=BinImage(fftshift(abs(fsum)),bin);
            % %             pspect0=BinImage(fftshift(abs(fmSum)),bin);
            %             subplot(2,3,2);
            % %             simg=Crop(pspect0,n/(bin*2)).^spectrumDisplayExp;
            %             imaga(trackDisDat.sp1.image);
            %             title('Before align: to 1/2 Nyquist');
            %
            %             subplot(2,3,3);
            %             simg=Crop(pspect,n/(bin*2)).^spectrumDisplayExp;
            %             imac(imscale(simg,255,spectrumDisplayLims));
            %             title('After align');
            %
            %             subplot(2,3,6);
            %             spect1d=Radial(pspect);
            %             spect1d(1)=NaN;
            %             plot((0:bin:n(1)/2-1)/(n(1)*pixA),sqrt(spect1d));
            %             ylabel('sqrt(spectrum)');
            %             xlabel('frequency, A^{-1}');
        end;
        %%  Write the output
        %             if numDefs>1
        %                 segName=[nameSegment num2str(segIndex)];
        %             else
        %                 segName='';
        %             end;
        outName=[mi.baseFilename pars.nameAli char(pars.nameSegment+segIndex-1)];
        if ~exist(mi.imagePath,'dir')
            mkdir(mi.imagePath);
        end;
        if pars.writeZTiff
            ztPars=struct;
            ztPars.snrRatio=200;
            ztPars.lfCutoff=.07;
            outputName=[mi.imagePath outName 'z.tif'];
            mdisp(pars.logs,['Writing ' outputName]);
            WriteZTiff(imgSum,pixA,outputName,ztPars);
            mi.imageFilenames{segIndex}=[outName '.mrc']; % sum filename is always stored as .mrc
        end;
        if pars.writeMRC
            mdisp(pars.logs,['Writing ' mi.imagePath outName '.mrc']);
            WriteMRC(imgSum,pixA,[mi.imagePath outName '.mrc']);
            mi.imageFilenames{segIndex}=[outName '.mrc']; % sum filename is always stored as .mrc
        end;
        WriteJpeg(BinImage(imgSum,2),[pars.dirJpeg outName '.jpg']);
        if pars.writeStack
            stackName=[mi.baseFilename pars.nameAstk char(pars.nameSegment+segIndex-1)];
            disp(['Writing the frame stack: ' stackName]);
            WriteMRC(rFrames,pixA,[mi.imagePath stackName '.mrc']);
        end;
    end;
    if displayOn && pars.doSaveFigures && pars.fig1>0
        set(gcf,'paperpositionmode','auto');
        jName=[pars.dirJpeg outName '-align.jpg'];
        print('-djpeg','-r200',jName);
    else
        save([pars.dirJpeg outName '-trackDisDat.mat'],'disDat');
    end;
end;  % for segIndex

end

function s=ShowFrameShifts(accumShiftsX, accumShiftsY, inds, iter, displayOn,s)
nim=size(accumShiftsX,1);
ymx=max([-accumShiftsX(:); accumShiftsY(:)]);
ymn=min([-accumShiftsX(:); accumShiftsY(:)]);
yspan=(ymx-ymn)*.1;

s.shiftX.linespec='.-';
s.shiftX.ploty=accumShiftsX(:,1:iter);
s.shiftX.axis=[0 nim+1 ymn-yspan ymx+yspan];
s.shiftX.yLabel='Original pixels';
s.shiftX.xLabel='Frame no.';
s.shiftX.title='X Shift';

s.shiftY=s.shiftX;
s.shiftY.ploty=accumShiftsY(:,1:iter);
s.shiftY.title='Y Shift';

if displayOn % Do the plotting directly
    subplot(234);  % show the cumulative translations
    %                 plot(-dsw*dxAccum(:,1:iter),'.-','markersize',5);
    plot(s.shiftX.ploty,'.-');
    axis(s.shiftX.axis);
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
    plot(s.shiftY.ploty,'.-');
    axis(s.shiftY.axis);
    title('Y Shift');
    set(gca,'xtick',1:nim,'xticklabel',xtext,...
        'xgrid','on','ygrid','on');
    xlabel('Frame no.');
    
    drawnow;
end;

end