% function TestEM12
% TestEM12
% Tests the reEMStep21 code

% Simulation parameters
% pixA=5;
% n0=64;
% angleSteps=[40 5 10];
% sigmaN=.1;

pixA=10;
n0=32;
angleSteps=[10 5 10];
% angleSteps=[10 10 20];
sigmaN=.003;
% sigmaN=.15;
% sigmaN=.1;  % n=64;
% sigmaN=.003;
nReps=[1];  % 10 images with each angle set
fscFlag=2;  % 0: compare with crystal; 1: normal FSC; 2: gold standard
nIters=4;
maxShift=3;
removeShifts=0;
useEM1=false;

symmetry=2;
s0=.001;
simAmp=.1;
simSigmaAmp=0;
simSigmaC=1;
simSigmaG=0;

imgAmp=.1;
sigmaC=1;
sigmaG=10;
dfMin=1;
dfRange=1;   % defocus from 1 to 2
nCTFs=10;
simbx0=0;     % image offset
simby0=0;
B0=50;     % B factor = B0 + B1*defocus
B1=50;
angleMins=[0 30 180];
simAngleMins=[0 30 180];
fc=pixA/60;  % make the starting volumes 60A resolution
trueParVals=[sigmaN^2 1 sigmaC^2 simby0 0];
parLabels={'noise var '
    'image amp '
    'click var '
    'b0        '
    '          '};

% fc=1;
nVols=numel(nReps);



% First, get the volumes.
origVols=arGetRefVolumes(pixA,n0,nVols);
refVols=zeros(size(origVols),'single');
for iVol=1:nVols
    refVols(:,:,:,iVol)=s0*SharpFilt(origVols(:,:,:,iVol),fc,fc/n0);
end;
% refVols=s0*origVols; %%%%%%%%%%
% refVols0=s0*origVols; %%%%%%%%%%

refVols=repmat(refVols,1,1,1,1,2);  % copy for twin reconstructions



% -------Make the fake data--------
si=struct;
si.pixA=pixA;

[simRi, simAngles]=reSetRefAngles(angleSteps,simAngleMins,1);  % include alphas
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

% si.ctfs=si.ctfs*0+1;%%%%%%%%%%----------------------------------------

si.miIndex=single(mod((0:nImgs-1)',nCTFs)+1);

clickX=simSigmaC*randn(nImgs,1)+simbx0;
clickY=simSigmaC*randn(nImgs,1)+simby0;
geomY=simSigmaG*randn(nImgs,1);

% Apply the ctfs and shifts to the images
imgs=zeros(size(imgsOrig),'single');
imgAmps=simAmp+simSigmaAmp*randn(nImgs,1);
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


moi.sigmaN=sigmaN;
moi.sigmaC=sigmaC;
moi.sigmaG=sigmaG;
moi.b0=0;
moi.pVols=ones(nVols,1)/nVols;
moi.imgAmps=imgAmp*ones(nImgs,1);
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

ri=struct;
ri.radius=floor(n0/2)-maxShift;
ri.nTrans=2*maxShift+1;
ri.softMask=fuzzymask(n0,2,ri.radius*.9,ri.radius*.2);

% Remove DC from images
invSoftMask=1-ri.softMask;
invSoftMask=invSoftMask/sum(invSoftMask(:)); % sums to 1
for i=1:size(imgs,3)
    img0=imgs(:,:,i);
    imgDC=invSoftMask(:)'*img0(:);
    imgs(:,:,i)=ri.softMask.*(img0-imgDC);
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
    [ri, refAngles]=reSetRefAngles(angleSteps,angleMins,false,ri);
    nRefs=size(refAngles,1);
    disp(['nAngles: ' num2str(ri.angleN) '  nRefs: ' num2str(nRefs)]);
    
    allProjs=zeros(n0,n0,nRefs,2,nVols,nTwins,'single');
    
    ri.radius=floor(n0/2)-maxShift;
    ri.nTrans=2*maxShift+1;

    ri.imgAmps=moi.imgAmps;
    ri.pVols=moi.pVols;
    ri.sigmaC=moi.sigmaC;
    ri.sigmaG=moi.sigmaG;
    ri.sigmaN=moi.sigmaN;
    
    if iter==2 && removeShifts  % remove the shifts
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
    
    
    if fscFlag==1  % conventional fsc
        refVols=repmat(mean(refVols,5),1,1,1,1,2); % force the two refVols to be identical.
    end;
    
        
    tic
    for iTwin=1:nTwins
        
        refs=reMakeTemplates(refVols(:,:,:,:,iTwin),refAngles);  % refs(x,y,iRef,iVol)
        figure(2);
        ImagicDisplay1(refs);
        
        ri.imgActive=twinFlags(:,iTwin);
        ri.accumClassMeans=false;
        % ------------------------ Do the EM step ---------------------------
        % ------------------------------------------------------------------
        disp(['EM step ' num2str(iTwin)]);
        if useEM1
            [classMeans, classNorms, ro1, pTrans1, pAlphas1, pRefs1]=reEMStep1(imgs,refs,si,ri);
        else
            [classMeans, classNorms, ro1, pTrans1, pAlphas1, pRefs1]=reEMStep22(imgs,refs,si,ri);
        end;
        % ------------------------------------------------------------------
        % ------------------------------------------------------------------
return        
        %%
        if iTwin==1
            roi=ro1;
            pTrans=pTrans1;
            pAlphas=pAlphas1;
            pRefs=pRefs1;
        else
            roi=AddStructFields(roi,ro1);
            pTrans=pTrans+pTrans1;
            pAlphas=pAlphas+pAlphas1;
            pRefs=pRefs+pRefs1;
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
%              allProjs(:,:,:,1,iVol,iTwin)=classMeans(:,:,:,iVol);
            allProjs(:,:,:,2,iVol,iTwin)=realNorms;
        end;
    end;
    toc

    %% ----------------Do the Fourier reconstruction----------------
    disp('Reconstructions');
    allFVols=reFourierInsertion(allProjs,refAngles);
%%    
    %     Normalization of volumes
    k=.000100;
    for iTwin=1:nTwins
        for iVol=1:nVols
            newVols(:,:,:,iVol,iTwin)=rsNormalizeReconstruction(allFVols(1,iVol,iTwin),allFVols(2,iVol,iTwin),k);
        end;
    end;
    toc
    
    %     figure(11);  % Show the 3D Fourier norms
    %     ShowSections2(allFVols(2).PadFT);
    %     subplot(3,3,1);
    %     title('Fourier norms');

    freqs=(0:n0/2-2)'/(n0*pixA);  % frequencies for fscs
    
    % Display volumes and fscs
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
        title('FSC');
        outputVols(:,:,:,iVol,iTwin,iter)=newVols(:,:,:,iVol,iTwin);
    end;
    
    %     for i=1:size(roi.ePars,2)
    %         str=sprintf('%s  %7.4f  %7.4f',parLabels{i},mean(roi.ePars(:,i)),trueParVals(i));
    %         disp(str);
    %     end;
    
    figure(4);
    SetGrayscale;
    imacs(mean(pTrans,3));
    
    % Assign new model values
    refVols=newVols;    
    moi.sigmaN=sqrt(mean(roi.varNs));
    moi.a=mean(roi.imgAmps);
    %     moi.sigmaC=sqrt(mean(roi.sigmaC));
    % sigmaG remains
    %     moi.b0=mean(roi.ePars(:,4));
    %     moi.pVols=mean(roi.pVols);
    moi.imgAmps=0*roi.imgAmps+mean(roi.imgAmps);
    
    moi
    %%
    figure(8);
    %     plot(roi.pVols(:,1));
    %     title('pVols');
    
    SetGrayscale;
    imacs(mean(pTrans,3));
    
    
    formatString='yyyymmdd-HHMMSS';
    outName=['/Users/fred/EMWork/Sims/TestEM' datestr(now,formatString) 'i' num2str(iter) '.mat'];
    disp(['Writing ' outName]);
    save(outName,'ri','moi','roi','newVols','pTrans','pAlphas','pRefs');
    %%
    % q=reshape(pAngles,ri.angleN(1),nRefs,nImgs);
    
    %     Show the angle probabilities
    figure(3);
    %     colormap jet
    SetGrayscale;
    subplot(1,1,1);
    %     iAlpha=ceil((size(pAngles,1)+1)/2);
    q=reshape(pRefs,nRefs*nVols,nImgs);
    imagesc(log(max(1e-8,q))/log(10));
    colorbar;
    title(['Iteration ' num2str(iter)]);
    
    figure(4);
    plot(roi.imgAmps);
    
    shifts=zeros(nImgs,2);
    if maxShift>0
        for i=1:nImgs
            [val, ix, iy]=max2di(pTrans(:,:,i));
            shifts(i,:)=[ix iy]-ceil((n0+1)/2);
        end;
    end;
    % plot(squeeze(sum(pAngles,2))');
end; % while iter
