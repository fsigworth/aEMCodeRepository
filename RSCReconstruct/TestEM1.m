% function TestEM1
% TestEM1

% Simulation parameters
% pixA=5;
% n0=64;
% angleSteps=[40 5 10];
% sigmaN=.1;

pixA=10;
n0=32;
angleSteps=[10 10 20];
sigmaN=.03;
nReps=[2 1];  % 10 images with each angle set
fscFlag=1;  % 0: compare with crystal; 1: normal FSC; 2: gold standard
nIters=1;

symmetry=2;
nVols=2;
s0=.001;
sigmaAmp=0;
sigmaC=2;
sigmaG=0;
dfMin=1;
dfRange=1;   % defocus from 1 to 2
nCTFs=10;
b0=2;     % image offset
B0=50;     % B factor = B0 + B1*defocus
B1=50;
angleMins=[-20 30 180];

fc=pixA/100;  % make the starting volumes 60A resolution
trueParVals=[sigmaN^2 1 sigmaC^2 b0 0];
parLabels={'noise var '
    'image amp '
    'click var '
    'b0        '
    '          '};

% fc=1;

% First, get the volumes.
origVols=arGetRefVolumes(pixA,n0,nVols);
refVols=zeros(size(origVols),'single');
for iVol=1:nVols
    refVols(:,:,:,iVol)=SharpFilt(origVols(:,:,:,iVol),fc,fc/n0);
end;
refVols=repmat(refVols,1,1,1,1,2);  % copy for twin reconstructions

% -------Make the fake data--------
si=struct;
si.pixA=pixA;

[simRi, simAngles]=reSetRefAngles(angleSteps,angleMins,1);  % include alphas
% If we want angles for each image, start with this:
% allSimAngles=[];
% nVolImgs=cumsum(nReps)*size(simAngles,1);
% for iVol=1:nVols
%     allSimAngles=[allSimAngles ; repmat(simAngles,nreps(iVol),1)];
% end;
%
nSimAngles=size(simAngles,1);
nImgs=0;
imgsOrig=zeros(n0,n0,nSimAngles*sum(nReps(1:nVols)),'single');
allSimAngles=zeros(nSimAngles*sum(nReps(1:nVols)),3);
for iVol=1:nVols
    nr=nReps(iVol);
    imgsTemp=reMakeTemplates(origVols(:,:,:,iVol),simAngles);
    nAngs=size(imgsTemp,3);
    for j=1:nAngs
        imgsOrig(:,:,nImgs+1:nImgs+nr)=repmat(imgsTemp(:,:,j),1,1,nr);
        allSimAngles(nImgs+1:nImgs+nr,:)=repmat(simAngles(j,:),nr,1);
        nImgs=nImgs+nr;
    end;
end;
imgsOrig=reshape(imgsOrig,n0,n0,nImgs);
%
% Make the ctfs and assign the indices
si.ctfs=zeros(n0,n0,nCTFs,'single');
for i=1:nCTFs
    defocus=dfMin+(i-1)/nCTFs*dfRange;  % step increase in def with index
    si.ctfs(:,:,i)=abs(CTF(n0,pixA,.025,defocus,2,B0+B1*defocus,.05));
end;
% si.ctfs=si.ctfs*0+1;%%%%%%%%%%
si.miIndex=single(mod((1:nImgs)',nCTFs)+1);

clickX=sigmaC*randn(nImgs,1);
clickY=sigmaC*randn(nImgs,1)+b0;
geomY=sigmaG*randn(nImgs,1);

% Apply the ctfs and shifts to the images
imgs=zeros(size(imgsOrig),'single');
imgAmps=1+sigmaAmp*randn(nImgs,1);
for i=1:nImgs
    p=FourierShift(n0,[clickX(i) clickY(i)]);
    imgs(:,:,i)=real(ifftn(fftn(imgsOrig(:,:,i)).*p.*...
        ifftshift(si.ctfs(:,:,si.miIndex(i)))))...
        *s0*imgAmps(i)+sigmaN*randn(n0,n0);
end;

si.rVesicle=single(200/pixA*ones(nImgs,1));  % vesicle radii are all the same
si.yClick=single(si.rVesicle.*sind(allSimAngles(:,2))+geomY-clickY);

fprintf('%d images\n',nImgs);
%%%
%  We now have imgs and the si struct initialized.
figure(1);
ImagicDisplay1(imgs);
drawnow;


moi.a=1;
moi.sigmaN=sigmaN;
moi.sigmaC=2;
moi.sigmaG=2;
moi.b0=0;
moi.pVols=ones(nVols,1)/nVols;
moi.imgAmps=imgAmps;
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
            %             classMeanSym=reSymmetrizeClasses(classMeans,ri,symmetry);
            allProjs(:,:,:,1,iVol,iTwin)=classMeans(:,:,:,iVol);
            allProjs(:,:,:,2,iVol,iTwin)=realNorms;
        end;
    end;
    toc
    %%     Do the Fourier reconstruction
    allFVols=reFourierInsertion(allProjs,refAngles);
    
    figure(11);
    ShowSections2(allFVols(2).PadFT);
    
    % %     Normalization of volumes
    k=.1;
    freqs=(0:n0/2-2)'/(n0*pixA);
    %%
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
    %%
    figure(8);
    plot(roi.pVols(:,1));
    title('pVols');
    
    formatString='yyyymmdd-HHMMSS';
    outName=['/Users/fred/EMWork/Sims/TestEM' datestr(now,formatString) 'i' num2str(iter) '.mat'];
    disp(['Writing ' outName]);
    save(outName,'ri','moi','roi','newVols','pTrans','pAngles');
    %%
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
