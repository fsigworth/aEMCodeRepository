function mi=k2FindJumpFcn(mi,pars,cPars)
%     function mi=k2FindJumpFcn(mi,pars,ctPars)
% Search for defocus jumps within a movie.
% The mi file is assumed to have the paths and movieName initialized.
% In this function the fields are set: mi.frameSets, mi.pixA and mi.dose
% Example of pars structure:
%     pars.defoci=[1.5 10];
%     pars.testSegs=[2 16; 30 inf];
%     pars.ds=4;
%     pars.ds1=4;
%     pars.fo=1/10;
%     pars.fi=1/200;
%     pars.doSaveFigures=1;
%     pars.showGraphics=0;
%     pars.dirJpeg='Jpeg/';
showPowerSpectra=0;

movieName=[mi.moviePath mi.movieFilename];
mdisp(pars.logs,['Reading ' movieName]);
[m,s]=ReadMovie(movieName);
pixA=mi.pixA;
[nx,ny,nim]=size(m);

testSegs=pars.testSegments;

% prune the test segments
nSegs=size(testSegs,1);
for i=2:nSegs
    if testSegs(i,1)>nim % a segment starts beyond the given segment start
        nSegs=i-1;
        testSegs=testSegs(1:nSegs,:);
        break;
    end;
end;
testSegs=min(nim,testSegs);
mi.frameSets=testSegs;  % default value
%%
mb8=BinImage(m(:),16);
mev=median(mb8);  % median of main part of images.
mev=max(mev,1);
if mi.cpe<=1  % do this only if we are handling raw counts
mdisp(pars.logs,'Removing outliers');
%     clear mb8;
%     Pick an outlier threshold
ithresh=ceil(4*mev);
while PoissonProb(ithresh,mev)>1e-6
    ithresh=ithresh+1;
end;
locs=find(m>ithresh);
if numel(locs)>0
    disp(['Clipping a total of ' num2str(numel(locs)) ' outliers to the value ' num2str(ithresh)]);
end;
m(locs)=ithresh;
end;
% Pad the image
disp('Padding');
n0=NextNiceNumber(max(nx,ny));
nw=n0/pars.ds;  % working image size
mc=zeros(n0,n0,nim,'single');
for i=1:nim  % crop individual images to avoid memory overuse.
    mc(:,:,i)=Crop(single(m(:,:,i)),n0,0,mev);  %??? 2*mev???
end;
if any(mi.imageSize)<n0
    mi.imageSize=[n0 n0];
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
    if sum(win1(:))<.99*numel(win1)  % significant masking to be done
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

% Test to see if any jump is to be found
nSegs=size(testSegs,1);
if nSegs>1 % we have to find the jump
    %%
    
    % find the tentative CTFs
    ctfs=zeros(nw,nw,nSegs);
    sqAccumCTF=zeros(nw,nw);
    disp('Estimating CTFs');
  
    if pars.showGraphics
    figure(1);
    end;

    cPars.showGraphics=pars.showGraphics;  % we set only this ctfOption, and kV.
    cPars.kV=mi.kV;
    cPars.minRes=pars.minRes;
    cPars.maxRes=pars.maxRes;
    
    for j=nSegs:-1:1
        msum=mean(mc(:,:,testSegs(j,1):testSegs(j,2) ),3);
        [ctp,~,ctfDisDat]=meFitCTF(msum,mi,pars,cPars,j,0);

        if j==nSegs
            ctPars=ctp;
        end;
        ctPars(j)=ctp;
        c=CTF(nw,pixA,ctp);
        ctfs(:,:,j)=c;
        sqAccumCTF=sqAccumCTF+c.^2;
        mdisp(pars.logs,[num2str(testSegs(j,:)) ': defocus ' num2str(ctp.defocus)]);
%         if ~pars.showGraphics
           ctfDatName=[pars.dirJpeg mi.baseFilename 'ju' char(96+j) '-ctfDisDat.mat'];
           save(ctfDatName,'ctfDisDat');
%         end;
    end;
    fs=RadiusNorm(nw);
    % sqAccumCTF=abs(ctfs(:,:,1).*ctfs(:,:,2)); % use the cross-term
    weights=sqrt(fs).*sqrt(sqAccumCTF);
    
    disp('FTs');
    fmc=fft2(BinImage(mc,pars.ds1));  % FT of each frame, after binning by ds1
    
    %% Make the Fourier mask and apply it.
    fmd=Cropo(fmc,nw,1);  % crop the whole stack
    ro=pixA*n0/pars.minRes(2);  % Outer mask radius
    ro=min(ro,n0*0.45);
    ri=pixA*n0/pars.maxRes(1);  % inner mask radius
    fmask=fuzzymask(nw,2,ro,ro/10)-fuzzymask(nw,2,ri,ri/10);
    fmdm=fmd;
    for i=1:nim
        fmdm(:,:,i)=fftshift(fmd(:,:,i)).*fmask;
    end;
    % fmdm is zero-centered.
    
    disp('Correlations');
    nSegs=size(testSegs,1);
    segSums=zeros(nw,nw,nSegs);
    
    corr=zeros(nim,nSegs);
    A=zeros(nSegs,nSegs);
    
    for j=1:nSegs
        i1=testSegs(j,1);
        i2=min(nim,testSegs(j,2));
        inum=i2-i1+1;
        %         partSums=zeros(nw,nw,i2,'single');
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
    for i=1:nSegs
        for j=1:nSegs
            if i~=j
                q=segSums(:,:,i);
                r=segSums(:,:,j).*weights;
                A(i,j)=real(q(:)'*r(:))/nw^2;
            end;
        end;
    end;
    
    for i=1:nim
        for j=1:nSegs
            if corr(i,j)==0  % we haven't computed it yet
                q=segSums(:,:,j);
                p=fmdm(:,:,i).*weights;
                corr(i,j)=real(q(:)'*p(:))/(nw^2);
            end;
        end;
    end;
    
%     Find the transition point pt
    amps=corr/A;
    dAmps=amps(:,2)-amps(:,1);
    diffAmps=diff(dAmps);
    [~,pt1]=max(diffAmps);  % sharpest increase is at transition
%     pt1
    pt1=max(pt1,testSegs(1,1)+1); % Can't be before the first test segment
    dAmps(1:pt1-1)=-1;  % 
    pt=find(dAmps>0,1);
    %     transitionPeak=ptJump+1
    %     transitionFrame=pt
    
    nbin=4;
    ndis=200;
    sp1s=zeros(n0/2/nbin,nim);
    imin=max(1,pt-4);
    imax=min(nim,pt+3);
    
    if pars.showGraphics && showPowerSpectra
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
    end;
    
    disDat=struct;
    me=mean(mc(:,:,testSegs(1,1):testSegs(1,2)),3);
    disDat.mean1.image=GaussFilt(BinImage(me,4),.1);
    disDat.mean1.title=mi.movieFilename;
    
    for i=1:10
        j=max(1,min(nim,pt-6+i));
        fn=['spect' num2str(i+1)];
        q=BinImage(mc(:,:,j),4);
        q=q-mean(q(:));
        n1=size(q,1);
        sp=fftshift(abs(fftn(q)).^2);
        sp=sp.*(1-fuzzymask(n1,2,n1/30,n1/100));
        disDat.(fn).image=BinImage(sp,4).^.5;
        disDat.(fn).axis='off';
        if i==1
            disDat.(fn).xlabel='Spectrum to 1/4 Nyquist';
        end;
        if j==pt
            disDat.(fn).title=['---' num2str(j) '---'];
        else
            disDat.(fn).title=num2str(j);
        end;
        
    end;
    
    disDat.corr.plotx=1:size(amps,1);
    disDat.corr.ploty=amps;
    disDat.corr.xLabel='Frame';
    disDat.corr.yLabel='Correlation';
%     plot(1:size(amps,1),amps,'.-','markersize',10);
%     hold on;
%     plot([pt],amps([pt],:),'rs');
%     hold off;
%     xlabel('Frame');
%     ylabel('Correlation');
    
    if pars.showGraphics %% pars.writeGraphics
        figure(3);
        DrawFigureFromData(disDat,3,4);
        set(gcf,'paperpositionmode','auto');
        outName=[pars.dirJpeg mi.baseFilename 'jump.jpg'];
        disp(['Writing ' outName]);
        print('-djpeg','-r200',outName);
    end;
    if pars.saveGraphics
        save([pars.dirJpeg mi.baseFilename '-jumpDisDat.mat'],'disDat');
    end;
    mi.frameSets(1,1)=testSegs(1,1);
    mi.frameSets(1,2)=pt-2;
    mi.frameSets(2,1)=pt+1;
    mi.frameSets(2,2)=testSegs(2,2);
else
    disp(['Movie has ' num2str(nim) ' frames, no jump assumed.']);
end; % if nSegs>1

mi.pixA=pixA;
mi.frameDose=borderVal/(mi.cpe*mi.pixA^2);
