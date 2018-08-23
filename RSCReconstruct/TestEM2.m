% function TestEM2
% Try reconstructing actual data


% [fname, pa]=uigetfile('*si.mat','Select an si file');
% [rootPath, stackPath]=ParsePath(pa);
% if isnumeric(rootPath)  % user clicked Cancel
%     return
% end;
% cd(rootPath);  % level above the stack directory
% disp(['Reading ' fname]);
% load([stackPath fname]);  % Losd the si structure
% %
% [pa, nm]=fileparts(fname);
% dataName=[stackPath nm(1:end-1) 'tack.mrc'];  % actual stack data
% disp(['Reading ' dataName]);
% [imgs0, imgPixA]=ReadEMFile(dataName);

%%
if abs(imgPixA-si.pixA)/si.pixA>.001
    error(['Discrepancy in pixel sizes: ' num2str([imgPixA si.pixA])]);
end;
nst=size(imgs0,1);  % size of stack images
nImgs=size(imgs0,3);

ds=1.75;
pixA=si.pixA*ds;  % working pixel size
n0=40;  % working image size

imgs=DownsampleGeneral(imgs0,n0,1/ds,2,1);  % downsample the stack

angleSteps=[40 5 10];
angleMins=[0 30 180];
fscFlag=1;  % 0: compare with crystal; 1: normal FSC; 2: gold standard
nIters=6;

symmetry=2;
nVols=1;
s0=median(si.sVesicle);

moi.sigmaN=std(imgs(:));
moi.sigmaC=2;
moi.sigmaG=0;
moi.a=1;
moi.b0=1;
moi.pVols=ones(nVols,1)/nVols;
moi.imgAmps=ones(nImgs,1);

fc=pixA/100;  % make the starting volumes 60A resolution
trueParVals=[moi.sigmaN^2 1 moi.sigmaC^2 moi.b0 0];
parLabels={'noise var '
    'image amp '
    'click var '
    'b0        '
    '          '};

% Get the starting volumes.
origVols=arGetRefVolumes(pixA,n0,nVols);
refVols=zeros(size(origVols),'single');
for iVol=1:nVols
    refVols(:,:,:,iVol)=SharpFilt(origVols(:,:,:,iVol),fc,fc/n0);
end;
refVols=repmat(refVols,1,1,1,1,2);  % copy for twin reconstructions


figure(1);
ImagicDisplay1(imgs);
drawnow;

%%
moi.a=1;
moi.sigmaN=std(imgs(:));
moi.sigmaC=2;
moi.sigmaG=2;
moi.b0=0;
moi.pVols=ones(nVols,1)/nVols;
moi.imgAmps=ones(nImgs,1);
s0=.01;

% refVols has already been assigned.

fsc=zeros(n0/2-1,nIters,nVols);
for iVol=1:nVols
    figure(4+iVol);
    ShowSections2(refVols(:,:,:,iVol));
    subplot(3,3,9);
    freqs=(0:n0/2-2)/(n0*pixA);
    fsc(:,1,iVol)=FSCorr(refVols(:,:,:,iVol),origVols(:,:,:,iVol));
    plot(freqs,fsc(:,:,iVol));
end;


%% iterate
nTwins=1+single(fscFlag>0);
nRefVols=1+single(fscFlag>1);

outputVols=zeros(n0,n0,n0,nVols,nTwins,nIters);
newVols=zeros(n0,n0,n0,nVols,nTwins);

iter=0;

% Write out the input data
formatString='yyyymmdd-HHMMSS';
outName=['/Users/fred/EMWork/Sims/TestEM' datestr(now,formatString) 'i' num2str(iter) '.mat'];
disp(['Writing ' outName]);
save(outName,'si','imgs','origVols','refVols');

%%

while iter<nIters
    iter=iter+1
    
    % --------Pick angles for the references----------
    [ri, refAngles]=reSetRefAngles(angleSteps,angleMins);
    disp(['nAngles: ' num2str(ri.angleN)]);
    nRefs=size(refAngles,1);
    
    allProjs=zeros(n0,n0,nRefs,2,nVols,nTwins,'single');
    
    ri.imgAmps=moi.imgAmps;
    ri.n=n0;
    
    if iter==2   % remove the shifts
        for i=1:nImgs
            p=FourierShift(n0,-shifts(i,:));
            imgs(:,:,i)=real(ifftn(fftn(imgs(:,:,i)).*p));
        end;
    end;
    
    %     Pick alternating images
    twinFlags=true(nImgs,1);
    if nTwins>1
        twinFlags(1:nTwins:end,1)=false;
        twinFlags(:,2)=~twinFlags(:,1);
    end;
    
    % Do the EM step
    emFlags=true([ri.angleN(1) nRefs nVols nImgs]);
    disp('reEMStep');
    
    if fscFlag==1  % conventional fsc
        refVols=repmat(mean(refVols,5),1,1,1,1,2); % force the two refVols to be identical.
    end;
    
    tic
    for iTwin=1:nTwins
        
        refs=s0*reMakeTemplates(refVols(:,:,:,:,iTwin),refAngles);  % refs(x,y,iRef,iVol)
        figure(2);
        ImagicDisplay1(refs);
        
        ri1=ri;
        ri1.imgAmps=ri.imgAmps.*twinFlags(:,iTwin);
        
        [classMeans, classNorms, ro1, pTrans1, pAngles1]=reEMStep(imgs,refs,emFlags,si,ri1,moi);
        
        %%
        if iTwin==1
            roi=ro1;
            pTrans=pTrans1;
            pAngles=pAngles1;
        else
            roi=AddStructFields(roi,ro1);
            pTrans=pTrans+pTrans1;
            pAngles=pAngles+pAngles1;
        end;
        
        
        %%
        for iVol=1:nVols
            realNorms=zeros(n0,n0,nRefs,'single');
            nCls=size(classNorms,3);
            for i=1:nCls
                realNorms(:,:,i)=fftshift(real(ifftn(ifftshift(classNorms(:,:,i,iVol)))));
            end;
            classMeanSym=reSymmetrizeClasses(classMeans,ri,symmetry);
            allProjs(:,:,:,1,iVol,iTwin)=classMeanSym(:,:,:,iVol);
            allProjs(:,:,:,2,iVol,iTwin)=realNorms;
        end;
    end;
    toc
    %     Do the Fourier reconstruction
    allFVols=reFourierInsertion(allProjs,refAngles);
    
    figure(11);
    ShowSections2(allFVols(2).PadFT);
    
    % %     Normalization of volumes
    k=.1;
    freqs=(0:n0/2-2)'/(n0*pixA);
    %
    for iTwin=1:nTwins
        for iVol=1:nVols
            newVols(:,:,:,iVol,iTwin)=rsNormalizeReconstruction(allFVols(1,iVol,iTwin),allFVols(2,iVol,iTwin),k)/s0;
        end;
    end;
    
    for iVol=1:nVols
        figure(4+iVol);
        ShowSections2(mean(newVols(:,:,:,iVol,:),5));
        subplot(3,3,9);
        if fscFlag>0
            fsc(:,iter+1,iVol)=FSCorr(newVols(:,:,:,iVol,1),newVols(:,:,:,iVol,2));
            target=0.143;
        else
            fsc(:,iter+1,iVol)=FSCorr(newVols(:,:,:,iVol),origVols(:,:,:,iVol));
            target=0.5;
        end;
        plot(freqs,[fsc(:,:,iVol) 0*freqs+target 0*freqs]);
        outputVols(:,:,:,iVol,iTwin,iter)=newVols(:,:,:,iVol,iTwin);
    end;
    
    for i=1:size(roi.ePars,2)
        str=sprintf('%s  %7.4f  %7.4f',parLabels{i},mean(roi.ePars(:,i)),trueParVals(i));
        disp(str);
    end;
    
    figure(4);
    SetGrayscale;
    imacs(mean(pTrans,3));
    
    % Assign new model values
    refVols=newVols;
    
    moi.sigmaN=sqrt(mean(roi.ePars(:,1)));
    moi.a=mean(roi.ePars(:,2));
    moi.sigmaC=sqrt(mean(roi.ePars(:,3)));
    % sigmaG remains
    moi.b0=mean(roi.ePars(:,4));
    moi.pVols=mean(roi.pVols);
    moi.imgAmps=1+0*roi.imgAmps;
    
    moi
    %
    figure(8);
    plot(roi.pVols(:,1));
    title('pVols');
    
    formatString='yyyymmdd-HHMMSS';
    outName=['/Users/fred/EMWork/Sims/TestEM' datestr(now,formatString) 'i' num2str(iter) '.mat'];
    disp(['Writing ' outName]);
    save(outName,'ri','moi','roi','newVols','pTrans','pAngles');
    %
    % q=reshape(pAngles,ri.angleN(1),nRefs,nImgs);
    
    %     Show the angle probabilities
    figure(3);
    colormap jet
    subplot(1,1,1);
    iAlpha=ceil((size(pAngles,1)+1)/2);
    q=sum(pAngles(:,:,1,:),1);
    imagesc(log(max(1e-8,squeeze(q)))/log(10));
    colorbar;
    title(['Iteration ' num2str(iter)]);
    
    figure(4);
    plot(roi.imgAmps);
    
    shifts=zeros(nImgs,2);
    for i=1:nImgs
        [val, ix, iy]=max2di(pTrans(:,:,i));
        shifts(i,:)=[ix iy]-ceil((n0+1)/2);
    end;
    % plot(squeeze(sum(pAngles,2))');
end; % while iter
