% k2FindDefocusJump.m
% Attempt to find the frame where the defocus jump occurs in a K2 movie.
% 
% Sets the following fields in the mi file:
% pixA
% frameDose
% frameSets

gPixADefault=1.247;

gFo=1/10;  % outer Fourier limit
gFi=1/200; % inner Fourier limit, A^-1
gds=4;
gds1=4;  % downsampling before fft; for speed, make it =ds

gTestSegs0=[2 16; 30 inf];
gDefoci=[1.5 10];

gSaveFigures=true;
% gMakeMiFiles=true;

if ~exist('gNames','var') || ~exist('gBatchProcessing','var') || ~gBatchProcessing
    [gNames, pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, ginfoPath]=ParsePath(pa);
    cd(rootPath);
    if ischar(gNames)
        gNames={gNames};
    end;
end;
%%

gNFiles=numel(gNames);

if gNFiles>1
    disp(['Working on ' num2str(gNFiles) ' files.']);
end;

gTransFrames=zeros(gNFiles,1);

for gIndex=1:gNFiles
    clear -regexp ^[^g] % clear everything that doesn't start with g
    infoName=[ginfoPath gNames{gIndex}];
    load(infoName);
    disp([num2str(gIndex) ':  ' gNames{gIndex} ': Reading ' mi.movieFilename]);
    movieName=[mi.moviePath mi.movieFilename];
    [m,s]=ReadMovie(movieName);
    pixA=s.pixA;
    if pixA==0  % file didn't contain the information
        pixA=gPixADefault;
        disp(['pixA set to default value ' num2str(pixA)]);
    end;

    [nx,ny,nim]=size(m);
    testSegs=min(nim,gTestSegs0);
    
    %%
    mb8=BinImage(m(:),16);
    [nx,ny,nim]=size(m);
    testSegs=min(nim,gTestSegs0);
    
    %%
    mb8=BinImage(m(:),16);
    mev=median(mb8);  % median of main part of images.
    clear mb8;
%     Pick an outlier threshold
    ithresh=ceil(4*mev);
    while PoissonProb(ithresh,mev)>1e-6
        ithresh=ithresh+1;
    end;
    locs=find(m>ithresh);
    if numel(locs)>0
        disp(['Clipping a total of ' num2str(numel(locs)) ' outliers to the value ' num2str(ithresh)]);
    end;
    mc(locs)=ithresh;
% Pad the image
    disp('Padding');
    n0=NextNiceNumber(max(nx,ny));
    nw=n0/gds;  % working image size
    mc=zeros(n0,n0,nim,'single');
    for i=1:nim  % crop individual images to avoid memory overuse.
        mc(:,:,i)=Crop(single(m(:,:,i)),n0,0,mev*2);  %??? 2*mev???
    end;
%%    
    disp('Masking')
    %   Mask any beam shadows
    beamEdgeA=300;  % no. of Angstroms in border
    msk=Crop(ones(nx,ny,'single'),n0);
%
    for iSeg=1:size(testSegs,1)
        i1=testSegs(iSeg,1);
        i2=testSegs(iSeg,2);
        me=mean(mc(:,:,i1:i2),3);  % look at summed image
        mev=mean(me(:));
        nBin=max(8,4*ceil(1/sqrt(mev)));   % at least 64 counts per bin
        mb=BinImage(me,nBin);  % smooth the image anyway: average bu 
        win1=mb>mean(mb(:))/2; % find regions below 1/2 of avg intensity.
        if sum(win1)<.99*numel(win1)  % significant masking to be done
%          disp('masking');
          msk=msk.*(GaussFiltDCT(ExpandImage(mb>mev/2,nBin),0.2*pixA/beamEdgeA)>.95);
        end;
        %%
    end;
    % we now have msk as the intersection of good regions in both test
    % segments.
    sumMasked=0;
    for i=1:nim
        m1=mc(:,:,i).*msk;
        mc(:,:,i)=m1;
        sumMasked=sumMasked+sum(m1(:));
    end;
    borderVal=sumMasked/(nim*sum(msk(:)));
    for i=1:nim
        m1=mc(:,:,i);
        m1(~msk)=borderVal;
        mc(:,:,i)=m1;
    end;
    
    %%
    
    % find the tentative CTFs
    nsegs=size(testSegs,1);
    ctfs=zeros(nw,nw,nsegs);
    sqAccumCTF=zeros(nw,nw);
    disp('Estimating CTFs');
    figure(1); SetGrayscale;
    for j=1:nsegs
        msum=mean(mc( :,:,testSegs(j,1):testSegs(j,2) ),3);
        [ctp, spectrum]=meFitCTF(msum,pixA,gDefoci(j),0);
        if j==1
            ctPars=ctp;
        else
            ctPars(j)=ctp;
        end;
        c=CTF(nw,pixA,ctp);
        ctfs(:,:,j)=c;
        sqAccumCTF=sqAccumCTF+c.^2;
        disp([num2str(testSegs(j,:)) ': defocus ' num2str(ctp.defocus)]);
    end;
    fs=RadiusNorm(nw);
    % sqAccumCTF=abs(ctfs(:,:,1).*ctfs(:,:,2)); % use the cross-term
    weights=sqrt(fs).*sqrt(sqAccumCTF);
    
    disp('FTs');
    fmc=fft2(BinImage(mc,gds1));  % FT of each frame, after binning by ds1
    
    %% Make the Fourier mask and apply it.
    fmd=Cropo(fmc,nw,1);  % crop the whole stack
    ro=gFo*pixA*n0;  % Outer mask radius
    ro=min(ro,n0*0.45);
    ri=gFi*pixA*n0;  % inner mask radius
    fmask=fuzzymask(nw,2,ro,ro/10)-fuzzymask(nw,2,ri,ri/10);
    fmdm=fmd;
    for i=1:nim
        fmdm(:,:,i)=fftshift(fmd(:,:,i)).*fmask;
    end;
    % fmdm is zero-centered.
    
    disp('Correlations');
    nsegs=size(testSegs,1);
    segSums=zeros(nw,nw,nsegs);
    
    corrs=zeros(nw,nw,nsegs,nsegs);
    corr=zeros(nim,nsegs);
    A=zeros(nsegs,nsegs);
    
    for j=1:nsegs
        i1=testSegs(j,1);
        i2=min(nim,testSegs(j,2));
        inum=i2-i1+1;
        partSums=zeros(nw,nw,i2,'single');
        segSums(:,:,j)=sum(fmdm(:,:,i1:i2),3);
        corrSum=0;
        for i=i1:i2
            q=fmdm(:,:,i);
            r=(segSums(:,:,j)-q).*weights;
            c=real(r(:)'*q(:))/((inum-1)*nw^2);
            corr(i,j)=c;
            corrSum=corrSum+c;
        end;
        A(j,j)=corrSum/inum;
        segSums(:,:,j)=segSums(:,:,j)/inum;
    end;
    %
    for i=1:nsegs
        for j=1:nsegs
            if i~=j
                q=segSums(:,:,i);
                r=segSums(:,:,j).*weights;
                A(i,j)=real(q(:)'*r(:))/nw^2;
            end;
        end;
    end;
    
    for i=1:nim
        for j=1:nsegs
            if corr(i,j)==0  % we haven't computed it yet
                q=segSums(:,:,j);
                p=fmdm(:,:,i).*weights;
                corr(i,j)=real(q(:)'*p(:))/(nw^2);
            end;
        end;
    end;
    
    
    amps=corr/A;
    dAmps=amps(:,2)-amps(:,1);
    diffAmps=diff(dAmps);
    pt=find(dAmps>0,1);
    [val,ptJump]=max(diffAmps);
    transitionPeak=ptJump+1
    transitionFrame=pt
    
    nbin=4;
    ndis=200;
    sp1s=zeros(n0/2/nbin,nim);
    imin=max(1,pt-4);
    imax=min(nim,pt+3);
    
    %     Make 1D power spectra
    for i=1:nim
        sp=RadialPowerSpectrum(BinImage(mc(:,:,i),nbin));
        sp1s(:,i)=GaussFilt(sp,.5);
    end;
    sp1s(:,imin)=mean(sp1s(:,2:imin),2);
    sp1s(:,imax)=mean(sp1s(:,imax:end),2);
    figure(2);
    semilogy(abs(sp1s(1:ndis,imin:imax)));
    legend(num2str((imin:imax)'));
    title([mi.movieFilename '  frame ' num2str(pt)],'interpreter','none')
    
    figure(3);
    nr=3;
    nc=4;
    SetGrayscale;
    subplot(nr,nc,1);
    me=mean(mc(:,:,testSegs(1,1):testSegs(1,2)),3);
    imacs(GaussFilt(BinImage(me,4),.1));
    subplot(nr,nc,2);
    title(mi.movieFilename,'interpreter','none');
    
    
    for i=1:10;
        j=max(1,min(nim,pt-6+i));
        subplot(nr,nc,i+1);
        q=BinImage(mc(:,:,j),8);
        q=q-mean(q(:));
        n1=size(q,1);
        sp=fftshift(abs(fftn(q)).^2);
        sp=sp.*(1-fuzzymask(n1,2,n1/30,n1/100));
        imacs(BinImage(sp,4).^.5);
        if j==pt
            title(['**' num2str(j) '**']);
        else
            title(j);
        end;
        
    end;

    subplot(nr,nc,12)
    plot(1:size(amps,1),amps,'.-','markersize',10);
    hold on;
    plot([pt],amps([pt],:),'rs');
    hold off;
    
    if gSaveFigures
        set(gcf,'paperpositionmode','auto');
        if ~exist('Jpeg','dir')
            mkdir('Jpeg');
        end;
        outName=['Jpeg/' mi.baseFilename 'jump.jpg'];
        disp(['Writing ' outName]);
        print('-djpeg','-r200',outName);
    end;
    mi.frameSets(1,1)=testSegs(1,1);
    mi.frameSets(1,2)=pt-2;
    mi.frameSets(2,1)=pt+1;
    mi.frameSets(2,2)=testSegs(2,2);
    frameSets=mi.frameSets
    mi.pixA=pixA;
    mi.frameDose=borderVal/(mi.cpe*mi.pixA^2);
    
%     mi.imageSize=n0;
    
    disp(['Updated: ' infoName]);
    disp(' ');
    save(infoName,'mi');

    gTransFrames(gIndex)=pt;
end; % for fIndex
% Show the final results

subplot(nr,nc,nr*nc-1);
hist(gTransFrames);
disp(' ');
disp('Summary:');
for gIndex=1:gNFiles
    disp([gNames{gIndex} '   ' num2str(gTransFrames(gIndex))]);
end;

clear m mc fm*  % get rid of the biggest variables.