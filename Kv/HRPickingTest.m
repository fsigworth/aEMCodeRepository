% HRPickingTest.m
% Locate particles on micrographs or in stacks
%

% % Fake data, start in this directory:
% cd ~/EMWork/Simulations/relion2/Fred2_HRPicking/
% starDir='Refine3D/job037/'
% starName='run_data.star'
% refDir=starDir;
% refName='run_class001.mrc'
% stackPrune=4;
% refCrop=96;
% bValue=0

% ----------
% W366F data start here:
% cd ~/EMWork/Yangyu/20210224_YW_sel/
starDir='Refine3D/job110/';
starName='run_data.star';
% To use the reconstructed map
% refDir='Postprocess/job171/';
% refName='postprocess.mrc';

% % To use the masked map
% refDir='';
% refName='postprocess_masked_2.55A.mrc';

% % To use our TM map
% refDir='HRPicking/';
% refName='tmMap.mrc';
% projsFilename='projsComp.mat';
%
% To use our alpha-subunit composite map
refDir='HRPicking/';
refName='compMap.mrc';
projsFilename='projsComp.mat';

bValue=80 % This gives the best results
stackPrune=0;
refCrop=216;
% -----------

computeNewProjs=1;
listStackFiles=0;
stackMode=1;
maxNumParticles=1e4; % actual set is 2678 particles.

skipStar=1;
if ~skipStar
    dataStarName=[starDir starName];
    disp(['Reading ' dataStarName '...']);
    [nms,dat]=ReadStarFile(dataStarName,1,maxNumParticles+200);
    op=dat{1};
    d=dat{2};
    disp('Decoding image names...');
    [d.imgNos,d.imgFiles]=rlDecodeImageName(d.rlnImageName);
    %%
    np0=numel(d.rlnImageName);
    [stackUNames,stackUPtrs,stackRPtrs]=unique(d.imgFiles);
    nUStacks=numel(stackUNames);
    [micUNames,micUPtrs,micRPtrs]=unique(d.rlnMicrographName);
    nUMics=numel(micUNames); % number of mic names = no. of CTFs
end;
%% Gather the whole stack of particles

stackFilename=stackUNames{1};
[stack1,sts]=ReadMRC(stackFilename);
n0=size(stack1,1); % size of the original stack
stackInds1=find(stackRPtrs==1);
np1=numel(stackInds1);
if listStackFiles
    disp(['Stack file: ' stackFilename '  ' num2str(np1) ' particles']);
end;
stack=stack1; % allocate the first bit
stackInds=stackInds1;
np=np1;

for iFile=2:nUStacks
    stackFilename=stackUNames{iFile};
    if ~exist(stackFilename,'file')
        disp(['File doesn''t exist: ' stackFilename])
        break;
    end;
    stack1=ReadMRC(stackFilename);
    stackInds1=find(stackRPtrs==iFile);
    np1=numel(stackInds1);
    if listStackFiles
        disp(['Stack file: ' stackFilename '  ' num2str(np1) ' particles']);
    end;
    stack(:,:,np+1:np+np1)=stack1(:,:,1:np1);
    np=np+np1;
    stackInds=[stackInds; stackInds1];
end;
disp(['done reading the stack. ' num2str(np) ' particles total.']);

%% prune the stack
npP=min(np,maxNumParticles);
if stackPrune>0 % Decimate the particle selection
    okInds=mod(stackInds,stackPrune)==1 & stackInds<=npP;
else
    okInds=true(npP,1);
end;
% Squash the inds and stack
% okInds(npP+1:end)=false; % take only npP particles before decimation
stackIndsP=stackInds(okInds); % pointers to original stack indices
stackP=stack(:,:,okInds);
npP=sum(okInds);

%% -----Compute all the projections-----
if computeNewProjs
    [tmRef,s]=ReadMRC([refDir refName]); % Get our 3D ref
    tmRef1=Crop(tmRef,refCrop);
    disp(['Making ' num2str(npP) ' projections']);
    [projs,angs,shifts]=rlMakeRelionProjections(tmRef1,d,stackIndsP,ct(1).pixA,100);
    npR=size(projs,3);
    disp(['Writing ' refDir projsFilename]);
    save([refDir projsFilename],'projs','angs','shifts','tmRef1');
else
    disp(['Reading ' refDir projsFilename]);
    load([refDir projsFilename]);
end;
%% Get a ctf for each micrograph
activeMicRPtrs=micRPtrs(stackIndsP);
[activeMics,activeMicLines]=unique(activeMicRPtrs);
disp('Making CTFs');
ct=rlStarLinesToCtf(nms,dat,activeMicLines);
nct=numel(ct);
ctfs=zeros(n0,n0,nct,'single'); % same size as particles
for i=1:nct
    ct(i).B=bValue;
    ctfs(:,:,i)=CTF(n0,ct(i));
end;

%% Operate on the projections with ctfs
npR=npP; % I guess it's just a name change.
cProjs=zeros(n0,n0,npR,'single');
for i=1:npR
    micInd=activeMicRPtrs(i);
    cProjs(:,:,i)= -real(ifftn(fftn(Crop(projs(:,:,i),n0)).*ifftshift(ctfs(:,:,micInd))));
end;
%% Show all the inner products
%
% nVec=n0^2;
% sVec=reshape(stackP(:,:,1:npP),nVec,npP);
% rVec=reshape(cProjs,nVec,npR);
% sqrs=diag(rVec'*rVec);
% rVecNorm=rVec./repmat(sqrs',nVec,1);
% cMat=(sVec'*rVecNorm);
% imacs(cMat);
%
% % imags(sVec'*sVec)
% return
%% Compute cross-correlations
ct0=ceil((n0+1)/2);
found=false(npR,1);
shiftErrs=zeros(npR,3);
figure(2);
displayOn=1;
listResults=1;

for i=1:npR
    q=stackP(:,:,i);
    pr=Crop(projs(:,:,i),n0);
    rr=cProjs(:,:,i);
    ccp=fftshift(real(ifftn(fftn(q).*conj(fftn(pr)))));
    
    cc=fftshift(real(ifftn(fftn(q).*conj(fftn(rr)))));
    [ccMax,sx,sy]=max2di(cc);
    rsx=sx-ct0;
    rsy=sy-ct0;
    
    if displayOn
        subplot(331);
        imags(GaussFilt(q,.2));
        title(i);
        
        subplot(332);
        imags(pr);
        
        subplot(333);
        imags(rr);
        
        subplot(335);
        imacs(ccp);
        
        
        subplot(336);
        imacs(cc);
        hold on;
        plot(sx,sy,'w+','linewidth',2)
        hold off;
        
        subplot(339);
        imacs(Crop(cc,16))
        hold on;
        plot(9,9,'wd','linewidth',2);
        plot(9+rsx,9+rsy,'w+','linewidth',2)
        hold off;
        colormap jet
        drawnow;
    end;
    
    
    shiftErrs(i,:)=[rsx rsy hypot(rsx,rsy)];
    if hypot(rsx,rsy)<2
        found(i)=true;
        marker=' ';
    else
        marker = 'x';
    end;
    if listResults || i==npR
        disp(['angs=  ' num2str(angs(i,:),'%8.1f') '  shifts=  ' num2str(shifts(i,:),'%8.1f')...
            '  ' num2str([rsx rsy],'%8.1f') '  correct:  ' sprintf('%d / %d %s', sum(found),i,marker)]);
    end;
end;
subplot(337);
plot([activeMicRPtrs found]);
subplot(338);
plot(shiftErrs(:,1),shiftErrs(:,2),'.');
return

%% Compute average projection spectrum
spSum=zeros(n0,n0,'single');
for i=1:npR
    spSum=spSum+abs(fftn(Crop(projs(:,:,i),n0)));
end;
%%
spMean=fftshift(spSum/npR);
spM1=Radial(spMean);
freqs=(1:n0/2)/(n0*s.pixA);
figure(3);
subplot(221);
imags(spMean);
subplot(222);
semilogy(freqs,spM1)
subplot(223);
plot(freqs,spM1.*(1:n0/2)'*2*pi);

spM3=RadialPowerSpectrum(tmRef);
subplot(224);
spM3d=spM3.*(1:n0/2)'*2*pi;
plot(freqs,spM3d);

%% Find signal as a function of defocus
% Alternate version, where we don't multply by diameter, giving an
% approximate SNR.
spM3d=spM3*200; %% don't multiply by r.
spM3d(1:8)=0;
figure(6)
defs=[.05:.01:1];
defs=[.08 .19 .31];
nDefs=numel(defs);
signals=zeros(nDefs);
ctPars=ct(1);
ctPars.B=80;
subplot(223);
plot(freqs,spM3d*20,'b-','linewidth',2);
hold on;
for i=1:nDefs
    ctPars.defocus=defs(i);
    c=ContrastTransfer(freqs,ctPars).^2;
    plot(freqs,c);
    signals(i)=c*spM3d;
end;
hold off;
% axis([0 inf 0 1])
legend('ref','.08','.19','.31')
% subplot(222);
% plot(defs,signals);

%% Look at the individual radial spectra of the projections
figure(5); clf;
sp1s=RadialPowerSpectrum(Crop(projs,n0,1),1);
%%
sp1sd=sp1s.*repmat(freqs',1,npP);
plot(freqs,sp1sd);
%% sort angles by tilt
[vals,inds]=sort(abs(angs(:,2)-90));
for i=1:4
    lower=floor((i-1)*npP/4+1);
    upper=floor(i*npP/4)
    subplot(2,2,i);
    plot(freqs,sp1sd(:,lower));
end;
%%
for i=1:npP
    plot(freqs,sp1sd(:,i));
    title(vals(i));
    drawnow;
end;

return



%% Try finding particles on micrographs

nMic=numel(micUNames);
for i=1
    if ~exist(micUNames{i},'file')
        continue;
    end;
    [mic,s]=ReadMRC(micUNames{i});
    micParticles=find(micRPtrs==i);
    np=numel(micParticles);
    if np<1
        continue;
    end;
    fMic=GaussFilt(mic,.05);
    [bX,bY,tX,tY]=MakeBoxDrawingVectors( ...
        [d.rlnCoordinateX(micParticles) d.rlnCoordinateY(micParticles)], ...
        op.rlnImageSize(1)/2,0.8);
    tStrings=cell(np,1);
    for j=1:np
        tStrings{j}=num2str(micParticles(j));
    end;
    figure(10);
    imags(fMic);
    hold on;
    plot(bX,bY,'color',[1 1 0]);
    text(tX,tY,tStrings,'color',[1 1 0], ...
        'HorizontalAlignment','left',  'VerticalAlignment','top');
    hold off;
    title(i);
    drawnow;
    
    %% try to find particles
    for j=1:np
        cp=cProjs(:,:,j);
        cpx=Crop(cp,[s.mx s.my]);
        
        
        
        
        
        
    end;
end;
