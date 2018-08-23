% k2DriftTracker
% Create a summed image from a K2 movie, read as an MRC file
% containing raw counts, using a gain reference if available.
% The stored aligned-summed images are padded to 3840 x 3840.
% If there are multiple segments, they are stored separately.
%
% Directory structure
%   The existing movie is assumed to reside in a directory like this:
% experiment/movie_frames/sq01_1/Sep24_18.52.48.mrc
% To use this program with flReadAllFilesInSubdirectories set, you would select
% the movie_frames directory to start the automatic scanning of files.
% With that flag unset, you would select individual .mrc files.
% 
%   alignment data is stored back into the same directory:
% experiment/movie_frames/sq01_1/Sep24_18.52.48-shiftsd1.mat
% experiment/movie_frames/sq01_1/Sep24_18.52.48-shiftsd1.mat
% 
%   The aligned segments are stored in a Micrograph directory 2 levels up
% experiment/Micrograph/sq01_1_001Sep24_18.52.48d1.mrc
% experiment/Micrograph/sq01_1_001Sep24_18.52.48d2.mrc
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

flReadAllFilesInSubdirectories=0;
namePost='';
nameAli='sumAli';
nameSegment='d';
flagRestoreImages=1;  % Compute the full-size aligned images
dsw=4;            % Downsampling of working images
numIters=4;
fig1Size=[1000 700];  % Size of main figure, in pixels.
flagSegments=[1 20; 23 inf];
% flagSegments=[2 11; 18 27];  % first, second good stretches

dirInfo='../Info/';
dirAli='../Micrograph/';  % where the output files are written
dirJpeg='../Micrograph/Jpeg/';
flagSaveFigures=1;

% Have the user select some movie files
if flReadAllFilesInSubdirectories
    rootPath=uigetdir('.','Select a directory containing directories of movie files');
    if isnumeric(rootPath)
        return
    end;
    cd(rootPath);
    d=dir;
    nd=0;
    pathName={};
    for i=3:numel(d)
        if d(i).isdir
        nd=nd+1;
        pathName{nd}=d(i).name;
        end;
    end;
    
    if nd>0  % directories found
        disp([num2str(nd) ' directories found']);
    else
        error('No enclosed directories found');
    end;
    
else
    [fname, pathName]=uigetfile('*.mrc','Select movie files','multiselect','on');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    if ~iscell(fname)
        fname={fname};
    end;
    pathName={pathName};
    [rootPath, dirMovie]=ParsePath(pathName{1});
    cd(rootPath);
    pathName{1}=dirMovie;
end;

%% Pick up a gain reference which will be used for all movies
[refName, refPath]=uigetfile('*.dm4','Select a gain reference');
if ischar(refName) % user selected something
    disp(['Loading the gain reference ' refPath refName]);
    gainRef=ReadEMFile([refPath refName]);
else
    gainRef=1;
    disp('No gain reference used.');
end;

if ~exist(dirAli,'dir')
    mkdir(dirAli);
end;
if ~exist(dirJpeg,'dir')
    mkdir(dirJpeg);
end;

% Loop over subd
for paIndex=1:numel(pathName)
    dirMovie=pathName{paIndex};
    disp(['Directory: ' dirMovie]);
    if flReadAllFilesInSubdirectories  % pick up all the files
        d=dir(dirMovie);
        numFiles=0;
        fname={};
        for i=1:numel(d)
            [~,nm,ex]=fileparts(d(i).name);
            if strcmpi(ex,'.mrc') % we take all .mrc files
                numFiles=numFiles+1;
                fname{numFiles}=d(i).name;
            end;
        end;
    end;
    
    %%
    numFiles=numel(fname);
    disp(['Working on ' num2str(numFiles) ' files']);
    %%
    tic
    for findex=1:numFiles
        numDefs=max(1,size(flagSegments,1));
        for segIndex=1:numDefs  % defocus stretch
            %         We seem to have trouble with a memory leak. Hence clear variables
            clearvars -except f* dir* dsw name* num* gainRef segIndex pa*
% Set up the main display window
            figure(1);
            SetGrayscale;
            sz=get(0,'screensize');  % Root monitor size.
            fig1Size=min(fig1Size,sz(3:4)*0.8);  % Don't allow the figure to be too big.
            pos=sz(3:4)/2-fig1Size/2;
            set(gcf,'outerposition',[pos fig1Size]);  % Put the figure in the center of the root monitor.
            [pnm, baseFilename]=fileparts(fname{findex});
            
            % Construct the output name
            serNumber=sprintf('%03d',findex);
            if dirMovie(end)=='/'
                dirMovie(end)=[];
            end;
            outBaseName=[dirMovie '_' serNumber baseFilename];
            
            movName=[AddSlash(dirMovie) fname{findex}];
            disp(['Reading ' movName]);
            %     Get the header information
            [mv s]=ReadMRC(movName,1,0);
            inds=flagSegments(segIndex,:);
            inds=min(s.nz,inds);
            nim=inds(2)-inds(1)+1;
            disp([' frames ' num2str(inds)]);
            pixA=s.pixA
            
            %     Try to minimize memory use by using singles mostly
            m0=zeros(s.ny,s.nx,nim,'single');
            msk1=true(s.ny,s.nx);
            
            mfMean=0;
            for i=1:nim
                m=single(rot90(ReadMRC(movName,i+inds(1)-1,1),3));
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
            m0=Crop((m0-mfMean).*repmat(msk0,1,1,nim),n,1);  % pad the entire stack.
            
            % Compute all the fts
            disp('Computing FTs');
            fm0=zeros([n nim],'single');
            for i=1:nim
                fm0(:,:,i)=fft2(double(m0(:,:,i)));
            end;
            %%
            clf;
            subplot(2,3,1);
            imac(imscale(BinImage(sum(m0,3),8),256,1e-5));
            axis off;
            title(baseFilename,'interpreter','none');
            
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
            nw=n/dsw;
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
            dxAccum=zeros(nim,numIters);  % accumulated shifts for plots
            dyAccum=zeros(nim,numIters);
            
            fShifts=complex(single(ones([nw nim])));  % Fourier shift arrays
            wISpec=zeros([nw nim]);  % Fourier-shifted spectrum models
            cciLocal=zeros(15,15,numIters);  % 1st ccs for figure
            
            %
            disp('Tracking');
            for iter=1:numIters
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
                
                ymx=max(-dsw*[dxAccum(:);dyAccum(:)]);
                ymn=min(-dsw*[dxAccum(:);dyAccum(:)]);
                yspan=(ymx-ymn)*.1;
                
                subplot(234);  % show the cumulative translations
%                 plot(-dsw*dxAccum(:,1:iter),'.-','markersize',5);
                plot(-dsw*dxAccum(:,1:iter),'.-');
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
                plot(-dsw*dyAccum(:,1:iter),'.-');
                axis([0 nim+1 ymn-yspan ymx+yspan]);
                title('Y Shift');
                set(gca,'xtick',1:size(dxAccum,1),'xticklabel',xtext,...
                    'xgrid','on','ygrid','on');
                xlabel('Frame no.');
                
                drawnow;
                
                if iter>1 && max(abs([xsh; ysh]))<0.1/dsw % shift is less than 0.1 original pixel
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
            
            if flagRestoreImages
                disp('Getting full-sized images');
                %%
                fsum=zeros(n);
                fshSum=zeros(n);
                for i=1:nim
                    fSh=FourierShift(n,-dsw*[dxTotal(i) dyTotal(i)]);
                    fshSum=fshSum+fSh;
                    fsum=fsum+fm0(:,:,i).*fSh;
                end;
                %%
                if numel(gainRef)>1
                    gainRefSq=Crop(gainRef.*msk1,n,0,1);  %%% use msk0 here???
                else
                    gainRefSq=msk1;
                end;
                imgSum=gainRefSq.*real(ifftn(fsum))+nim*mfMean;  % restore the mean
                %
                spectrumDisplayExp=.1;
                spectrumDisplayLims=[0 5e-4];
                bin=16;
                figure(1);
                subplot(2,3,1);
                imac(imscale(BinImage(imgSum,4),256,1e-5));
                axis off;
                title(baseFilename,'interpreter','none');
                
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
                %%
                % Write the output
                if numDefs>1
                    segName=[nameSegment num2str(segIndex)];
                else
                    segName='';
                end;
                outName=[outBaseName nameAli segName];
                disp(['Writing ' dirAli outName '.mrc']);
                WriteMRC(imgSum,pixA,[dirAli outName '.mrc']);
                WriteJpeg(rot90(imgSum),[dirJpeg outName '.jpg']);
                disp('done.');
                disp('----');
                if flagSaveFigures
                    set(gcf,'paperpositionmode','auto');
                    jName=[dirJpeg outName '-align.jpg'];
                    print('-djpeg','-r200',jName);
                end;
                shiftX=dxTotal*dsw;
                shiftY=dyTotal*dsw;
                save([AddSlash(dirMovie) baseFilename '-shifts' segName '.mat'],'shiftX','shiftY');  % Save the working image set
            end;
            %%
            %
        end;
    end;
end;
toc;
