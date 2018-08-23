% function k2DriftTrackFcnQuick(movieName,pars)

% This is the quick version using GPU when possible, with the goal of
% showing the power spectrum only for the aligned movie.

% % test code
pars.useParfor=1;
pars.useGPU=1;
pars.segment=[1 20];
pars.dsw=2;
pars.ds0=2;
pars.pixA=1.2;
% cd('/Users/fred/EMWork/Hideki/140523Louise/movie_frames/sq02_2')
movieName='Jun10_23.26.15.tif';
% % end of test code


ds0=pars.ds0;  % binning of raw image
dsw=pars.dsw;  % further downsampling for tracking
pixA=pars.pixA*ds0;
numIters=2;

figure(1);
SetGrayscale;

tic
disp(['Reading ' movieName ]);
%     Get the header information
[~, s]=ReadMovie(movieName,1,0);
inds=pars.segment;
inds=min(s.nz,inds);
nim=inds(2)-inds(1)+1;
disp([' frames ' num2str(inds)]);

%     Try to minimize memory use by using singles mostly
m0=zeros(s.ny/ds0,s.nx/ds0,nim,'single');  % downsample raw images by ds0
msk1=true(s.ny/ds0,s.nx/ds0);
n=[1 1]*NextNiceNumber(max(s.nx,s.ny)/ds0);

mfMean=0;
m1=zeros([n nim],'single');
if pars.useParfor % reads from my SSD 3x faster in parallel, 2.5 instead of 7.5 (17 frames)
    parfor i=1:nim
        m=rot90(single(BinImage(ReadMovie(movieName,i+inds(1)-1,1),ds0)));
        %         m=single(rot90(ReadMovie(movName,i+inds(1)-1,1),3));
        %                 in case there are ADC failures, mask all vertical lines that are zeros
        msk0=repmat(any(m,2),1,s.nx/ds0);
        msk1=msk1&msk0;
        m=m.*msk0;
        pixMean=sum(m(:))/sum(msk0(:));
        m1(:,:,i)=Crop((m.*msk0)-pixMean,n);  % pad the entire stack.
        mfMean=mfMean+pixMean;
    end;
else
    for i=1:nim
        m=rot90(single(BinImage(ReadMovie(movieName,i+inds(1)-1,1),ds0)));
        %         m=single(rot90(ReadMovie(movName,i+inds(1)-1,1),3));
        %                 in case there are ADC failures, mask all vertical lines that are zeros
        msk0=repmat(any(m,2),1,s.nx/ds0);
        msk1=msk1&msk0;
        m=m.*msk0;
        pixMean=sum(m(:))/sum(msk0(:));
        m1(:,:,i)=Crop((m.*msk0)-pixMean,n);  % pad the entire stack.
        mfMean=mfMean+pixMean;
    end;
end;
toc

disp('Computing FTs');
if pars.useGPU
    %     % Compute all the fts on gpu
    fm0=zeros([n nim],'single');
    for i=1:nim % do it one by one, othwerwise not enough gpu memory...
        gm1=gpuArray((m1(:,:,i)));
        fm0(:,:,i)=gather(fftn(gm1));
    end;
else
    % Compute all the fts
    fm0=zeros([n nim],'single');
    for i=1:nim
        gm1=((m1(:,:,i)));
        fm0(:,:,i)=(fftn(gm1));
    end;
end;
toc
%%
subplot(1,2,1);
imac(imscale(BinImage(sum(m1,3),4),256,1e-5)); % Show the raw sum
title(movieName,'interpreter','none');
drawnow;

% Compute all the ccfs
fmSum=sum(fm0,3);
ncc=256;

ccs =zeros(ncc,ncc,nim,'single');
ctr=ncc/2+1;
disp('Computing overall CCF');
sumCs=zeros(n,'single');
for i=1:nim
    cs=fm0(:,:,i).*conj((fmSum-fm0(:,:,i)));
    sumCs=sumCs+cs;
end;
sumCs(1,1:120/ds0:n(2))=0;  % Remove the spikes
ccMean=Crop(fftshift(real(ifftn(sumCs))),ncc)/nim;

% Estimate the artifact amplitude

us=-.21;  % undershoot values !! should be modified for binning.
us2=-.08;
mul=[1 1 1 1 1; 1 0 0 0 1; 1 1 0 1 1; 1 0 0 0 1; 1 1 1 1 1]';
pmul=Crop(mul,ncc);
frac=[us us2 us; 0 1 0; us us2 us]';
pfrac=Crop(frac,ncc);
%     Estimate the value of the peak artifact
pkv=ccMean(ctr,ctr)-(pmul(:)'*ccMean(:))/sum(mul(:));
ccCorr=pkv*pfrac;
%%

% Downsample the images for alignment search
disp('Downsampling');
nw=n/dsw;
nwx=nw(1);
nwy=nw(2);
dsWindow=ifftshift(fuzzymask(nw,2,0.95*nw(1),.1*nw(1)));
wfCorr=dsWindow.*Cropo(fftn(ifftshift(Crop(ccCorr,n))),nw);

if pars.useGPU
    gwF0=zeros([nw nim],'single','gpuArray');
    gwI0=zeros([nw nim],'single','gpuArray');
    gdsWindow=gpuArray(single(dsWindow));
    for i=1:nim
        gwF0(:,:,i)=gdsWindow.*gpuArray(Cropo(fm0(:,:,i),nw));  % truncated ft
        gwI0(:,:,i)=real(ifftn(gwF0(:,:,i)));
    end;
    wF0=gather(gwF0);
    wI0=gather(gwI0);
    
else
    wF0=zeros([nw nim],'single');
    wI0=zeros([nw nim],'single');
    for i=1:nim
        wF0(:,:,i)=dsWindow.*Cropo(fm0(:,:,i),nw);  % truncated ft
        wI0(:,:,i)=real(ifftn(wF0(:,:,i)));
    end;
end;
%%  Initialize tracking
wI=wI0;
wF=wF0;
xsh=zeros(nim,1);
ysh=zeros(nim,1);
dxTotal=single(zeros(nim,1));
dyTotal=single(zeros(nim,1));
dxAccum=zeros(nim,numIters);  % accumulated shifts for plots
dyAccum=zeros(nim,numIters);

fShifts=complex(single(ones([nw nim])));  % Fourier shift arrays
wISpec=zeros([nw nim]);  % Fourier-shifted spectrum models
toc
%
disp('Tracking');
for iter=1:numIters
    disp(['Iteration ' num2str(iter)]);
    
    % Construct the independent means
    meanInds=1:nim;  % indices of images to sum for means.
    wFMean=zeros([nw nim],'single');
    for i=1:nim
        meanIndsX=meanInds(meanInds~=i);  % Delete one entry
        %             wIMean(:,:,i)=mean(wI(:,:,meanIndsX),3);
        wFMean(:,:,i)=mean(wF(:,:,meanIndsX),3);
        wISpec(:,:,i)=wfCorr.*fShifts(:,:,i).*conj(mean(fShifts(:,:,meanIndsX),3));
    end;
    toc
    %% Compute the weighting filter
    disp('Weighting filter');
    [H1, rs]=k2GetWeightingFilter5(wI,wFMean,wISpec,5,0);
    toc
    
    %% Compute the shifts from the cross-correlation peaks
    ccSum=zeros(nw);
    for i=1:nim
        crossSpec=H1.*wF(:,:,i).*conj(wFMean(:,:,i));
        cc=fftshift(real(ifftn(crossSpec-H1.*wISpec(:,:,i))));
        ccSum=ccSum+cc;
        [mxv, ix, iy]=max2di(cc);
        xsh(i)=ix-nwx/2-1;
        ysh(i)=iy-nwy/2-1;
    end;
    % Accumulate the shifts
    dxTotal=dxTotal+xsh;  % used for image correction
    dxAccum(:,iter)=dxTotal;     % 2d array for plotting
    dyTotal=dyTotal+ysh;
    dyAccum(:,iter)=dyTotal;
    %  Display the progress.
    ymx=max(-dsw*[dxAccum(:);dyAccum(:)]);
    ymn=min(-dsw*[dxAccum(:);dyAccum(:)]);
    yspan=(ymx-ymn)*.1;
    
    subplot(2,4,3);  % show the cumulative translations
    %                 plot(-dsw*dxAccum(:,1:iter),'.-','markersize',5);
    plot(-dsw*dxAccum(:,1:iter),'.-');
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
    set(gca,'xtick',1:size(dxAccum,1),'xticklabel',xtext,...
        'xgrid','on','ygrid','on');
    cellTxt=cell(iter,1);
    for j=1:iter
        cellTxt{j}=num2str(j);
    end;
    legend(cellTxt,'location','northeast');
    
    subplot(2,4,4);
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
    toc
    %% Shift the working images, and update the Fourier shifts
    disp('shifting')
    if pars.useGPU
        for i=1:nim
            gSh=gFourierShift(nw,-[dxTotal(i) dyTotal(i)]);
            %             gF0=gpuArray(wF0(:,:,i));
            gF=gpuArray(wF0(:,:,i)).*gSh;
            %             gF=gF0.*gSh;
            wF(:,:,i)=gather(gF);
            wI(:,:,i)=gather(real(ifftn(gF)));
            fShifts(:,:,i)=gather(gSh);
        end;
    else
        for i=1:nim
            fSh=FourierShift(nw,-[dxTotal(i) dyTotal(i)]);
            wF(:,:,i)=wF0(:,:,i).*fSh;
            wI(:,:,i)=real(ifftn(wF(:,:,i)));
            fShifts(:,:,i)=fSh;
        end;
    end;
    toc
end; % for iter

%%  Get the full-size spectrum

disp('Shifting whole image');

if pars.useGPU
    gsum=zeros(n,'single','gpuArray');
    for i=1:nim
        gSh=gFourierShift(n,-dsw*[dxTotal(i) dyTotal(i)]);
        gsum=gsum+gpuArray(fm0(:,:,i)).*gSh;
    end;
    fsum=gather(gsum);
else
    fsum=zeros(n);
    for i=1:nim
        fSh=FourierShift(n,-dsw*[dxTotal(i) dyTotal(i)]);
        fsum=fsum+fm0(:,:,i).*fSh;
    end;
end;
sp2=abs(fftshift(fsum)).^2;

nx=n(1);
xs=(-nx/2:nx/2-1)./(nx*pixA);
subplot(1,2,1);
imac(xs,xs,imscale(GaussFilt(sp2,.1).^.3,256,[0 2e-3]));
subplot(2,2,4);
s=Radial(sp2);
pt=round(nx/80);  % truncate the spectrum below this point
s(1:pt)=s(pt);
semilogy(xs(nx/2+1:nx),s);
toc
