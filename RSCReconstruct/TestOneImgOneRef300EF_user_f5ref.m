% TestOneImgOneRef300EF_user_f5.m
% derived from TestOneImgOneRef300EF_user.m
% by default run split reconstruction of a dataset
% load tsi.mat and read tstack.mrc to start
% frequently changed parameters are listed before stack input session
% can compute gold standard fsc with split data or traditional fsc with the
% whole dataset
% save:  
% (a) Maximum likelihood results 
% (b) last volume with membrane density added, fsc
% (c) unmask volume
% (d) sum of two volumes from split dataset, membrane added, gold standard
% fsc
% when compute the split dataset, wait until the split reconstruction
% completed, and run FSCgold_user
% Yi, Oct 29, 2013
%% parameter setup

pixA=5;

%%  load the intial model
mapName='/Volumes/TetraData/Structures/AMPAR/3KG2mapsub2.9A.mat';
load(mapName);
pixA0=2.9;  % model pixel size

nm=size(map,1);
map0=circshift(map,[0 0 3]);  % shift to center it.
ds=pixA/pixA0;

n=64;
map1=Crop(DownsampleGeneral(map0,round(nm/ds)),n); % resampled -> map0
ShowSections2(map1);

ri=reSetRefAngles([5 5 10],[10 30 180]);



return


display('niter dAlpha peakAlpha dBeta dGamma');
display([niter dAlpha peakAlpha dBeta dGamma]);
% parameter_input=input('Input a 1*5 vector in [] if you want to change them, otherwise press Enter: ');
% if parameter_input
%     [niter dAlpha peakAlpha dBeta dGamma]=parameter_input;
% end;
%refBetas=10:10:170;    % angle step  avoid beta=0 and 180 for b0
%refBetas=30:10:150;
%refGammas=0:15:165;
 refBetas=30:dBeta:150;
 refGammas=0:dGamma:170;
%refBetas=20:20:160;
%refGammas=0:30:165;
thre_gm=0.02; % amplitude threshold to generate mask
FlatModel=0;
dsi=1; % image downsample rate

%%  load image stack
if Fake
    n=40;
    disp('This is a simulation');
else
    disp('Operated on experimental images')
    %load sq03_10027tsi.mat;% 3124 particles with merged ctf
    %imgs0=ReadMRC('sq03_10027stall.mrc');
    %load ../RSC-ML2/sq03_10027w010e1tsi.mat;% 3124 particles with second exposure
    %imgs0=ReadMRC('../RSC-ML2/sq03_10027w010e1stall.mrc');
    [stackname, pa]=uigetfile('*stack.mrc','Select datastack');
    imgs0=ReadMRC([pa stackname]);
    nImgs0=size(imgs0,3);
    mIndices0=whichhalf:particle_split:nImgs0;
    [infoname, pa_info]=uigetfile('*tsi.mat','Select tsi info file');
    load([pa_info infoname]);
    mIndices0=whichhalf:particle_split:nImgs0;
    [indexname, pa_ind]=uigetfile('*Index.mat','Select index file if you need');
    if indexname
        load([pa_ind indexname]);    %load *Index.mat  % select certain particles 
        %nImgs0=numel(mIndicesInFredsStack);
        %mIndices0=whichhalf:particle_split:nImgs0;
        mIndices=mIndicesInFredsStack(mIndices0);
    else
        mIndices=mIndices0;
    end;    
    
    
    %load sq02_10177tsi.mat;
    %imgs0=ReadMRC('sq02_10177tstack.mrc');
    %load k1Index.mat; %load mIndicesInFredsStack
    %imgs0=ReadMRC('k1stack.mrc');  % this is katrine's stack, nImgs doesn't match
   
%    load sq02_1_2020Oct05_22.05.46as1tsi.mat; %************************  input file
%    imgs0=ReadMRC('sq02_1_2020Oct05_22.05.46as1tstack.mrc'); %*********  input file
    


    nImgs=numel(mIndices);
    pixA=si.pixA*dsi;
    y0s=si.yClick(mIndices)/dsi;
    vesRs=si.rVesicle(mIndices)/dsi;
    sVesicle=si.sVesicle(mIndices);
    ctfIndex=si.miIndex(mIndices);   
       
    %nImgs=numel(y0s); 
    n0=size(imgs0,1); 
    n=n0/dsi;
    ctfs_c=Crop(si.ctfs,n,1);
    nctfs=size(si.ctfs,3);
    if n==n0
        imgs=imgs0(:,:,whichhalf:particle_split:nImgs0);
    else
        display('downsample the images')
        imgs=Downsample(imgs0(:,:,whichhalf:particle_split:nImgs0),n,1);
    end;
    ImgMean=mean(imgs(:));
    if n0<(320-1)/si.pixA
        display('padding the images');
        n=round(320/pixA);
        imgs=Crop(imgs,n,1,ImgMean);
        ctfs=zeros(n,n,nctfs);
        for ictf=1:nctfs
            ctfs(:,:,ictf)=DownsampleGeneral(ctfs_c(:,:,ictf),n,n/(n0/dsi));
        end;
    else
        ctfs=ctfs_c;
    end;  
    % msk=fuzzymask(n,2,0.45*n,.1*n);  too tight for efd
    msk=fuzzymask(n,2,0.45*n,0.05*n);

    
    for i=1:nImgs  %katrine
        imgs(:,:,i)=msk.*(imgs(:,:,i)-ImgMean);
    end;

   
%     ntemp=round(n0*pixA0/pixA); % for downsampling
%     tempimg=zeros(ntemp,ntemp,nImgs);
%     for i=1:nImgs
%         tempimg(:,:,i)=DownsampleGeneral(imgs0(:,:,i),ntemp);
%     end;
%     imgs=Crop(tempimg,n,1);
    %imgs=imgs(n0/2-n/2+1:n0/2+n/2,n0/2-n/2+1:n0/2+n/2,:);
    
    
    for ictf=1:nctfs   
        ctfs(:,:,ictf)=ifftshift(ctfs(:,:,ictf));
    end;
    % exclude particles
    d0=160;  % angstroms min freqency
    d1=80; % angstroms
    sps=RadialPowerSpectrum(imgs);
        fmin=floor(pixA*n/d0);
        fmax=ceil(pixA*n/d1);
        lfv=(fmin:fmax)*sps(fmin:fmax,:);
        sigmap=std(lfv);
        threshSD=2.5;
    exc=lfv>threshSD*sigmap;
    sum(exc)
    %ctfIndex(exc)=0;
end;

%% initialize values and parameter


%FOR SIMULATION
ds=2;
a=.01;
vesR=100/ds;
nrep=6;


dfmin=2;
dfspan=1;

% a=.1;
sigmaN=1;
sigmaR=5;
sigmaC=2;
sigmaG=2;
b0=0;

npars=5;
theoret=[sigmaN^2 a sigmaC^2 b0 0];
%%%%%%%%%%%%%%%
aS(2)=0.4;  %initial guess of parameters
aS(1)=0.8;
%sigmaR=5;
aS(3)=8;
%sigmaG=2;
aS(4)=2;
%aS(4)=0.5;
aS=theoret;  % aS(5)  never used, could go into simulation
%%  load the intial model
%load /Users/fred/matlabWork/RSC/3KG2map5.8A.mat
%load /home/yc323/RSC-ML2/3KG2map5.8A.mat
%load /Users/yi/Installers/EMData/RSC201304/RSC/3KG2map5.8A.mat
load ../RSC2/3KG2RotMap2.9A.mat
pixA0=2.9;  % model pixel size
nm=size(map,1);
dm=3;
map(:,:,1+dm:nm)=map(:,:,1:nm-dm);  %SHIFT THE MAP by 7 to center 
map(:,:,1:dm)=0;
ds=pixA/(pixA0); %6.97 3.486is the pixel size for the stack
map0=Crop(DownsampleGeneral(map,round(nm/ds)),n); % resampled -> map0
mapx=map0; % intialize the size
%pixA=pixA0*ds;
mem_ptcl=50;  % in angstrom 180/2-40
membrane=round(n/2-mem_ptcl/pixA); 
%membrane=15;
%mapx(:,:,1+4:n)=mapx(:,:,1:n-4);
%mapx(:,:,1:n-3)=mapx(:,:,1+3:n);
 
fc=1/(60/pixA); % frequency threshold to generate volmsk and vol
if FlatModel % whether to flatten map0 above membrane region for mapx and Volmsk
    map_w=map0;
    for i=membrane:n
        map_w(:,:,i)=DownsampleGeneral(map0(:,:,i),n,1.2); %1.2 IS widened
    end;
   % map_wi=map_w;
    for i=1:n
        tempmap=DownsampleGeneral(squeeze(map_w(i,:,:)),n,1/1.2);
        map_w(i,:,membrane:n)=tempmap(:,membrane:n);
    end;

    mapx(:,:,1:membrane)=map_w(:,:,1:membrane)*2/3;  % membrane subtraction from mapw->mapx
    Volmsk=(GaussFilt(map_w,fc/2)>thre_gm*std(map_w(:)));
else 
    mapx(:,:,1:membrane)=map0(:,:,1:membrane)*2/3; % membrane subtracted from normal
    Volmsk=(GaussFilt(map0,fc/2)>thre_gm*std(map0(:)));
end;


Volmsk=GaussFilt(Volmsk,0.06);  % make it softer
vol=SharpFilt(mapx,fc,0); % the initial volume
RefScale=sum(vol,3);
maxRefScale=max(RefScale(:));
figure(2);
SetGrayscale;
ShowSections(mapx-Volmsk.*mapx);
drawnow;
%% prepare to make references
angs=[];
k=0;    

nRefBetas=numel(refBetas);
nRefGammas=numel(refGammas);
for beta=refBetas
    for gamma=refGammas;  % There are 6 gammas
        k=k+1;
        angs(k,:)=[0 beta gamma];
    end;
end;
%refMsks=rsMakeTemplates(angs,Volmsk);  this is not bia-free
%refs=GaussFilt(rsMakeTemplates(angs,mapx),0.05,1);
% refMsks0=refs>0.01*std(refs(:));
% refMsks=GaussFilt(refMsks0,0.05,1);

nRefs=k;
allBetas=angs(:,2);
betaRange=[min(allBetas,0) allBetas+5];
priorPKs=mean(sind(betaRange),2).*(betaRange(:,2)-betaRange(:,1))*pi^2/180^2;
%nRefBetas, % add  along each beta, to evaulate alphas fro all the images 

%refs=rsMakeTemplates(angs,mapx);

% figure(3);
% ImagicDisplay1(refs,2);

%% Create the fake images or modify parameter guessing
if Fake
imBetas=20:160;
%imBetas=[ 20 40 60 80];
imGammas=[0 30 60 90 120 150];
[rBetas rGammas]=ndgrid(imBetas,imGammas);
imAngs=repmat([rBetas(:)*0 rBetas(:) rGammas(:)],1,nrep);
nImgs=numel(rBetas)*nrep;

imAngs=reshape(imAngs',3,nImgs)';

rocks=sigmaR*randn(nImgs,2);
    rocks(:,2)=0;  % zero beta rocking
bobs=zeros(nImgs,1);
clicks=sigmaC*randn(nImgs,2);
clicks(:,2)=clicks(:,2)-b0;
[imgs0 y0s]=rsMakeFakeImages(mapx,0,vesR,imAngs,rocks,bobs,clicks);
ctfs=single(zeros(n,n,nImgs));

% operate with ctfs
for i=1:nImgs
    defocus=dfmin+rand*dfspan;
    ctfs(:,:,i)=ifftshift(abs(CTF(n,pixA,.025,defocus,2,100+100*defocus,.05)));
end;
ctfIndex=1:nImgs;
ctfIndex(1:2:5)=0;
imgsf=real(ifft2(fft2(imgs0).*ctfs));
vesRs=vesR*ones(nImgs,1);
msk=fuzzymask(n,2,0.45*n,.1*n);
imgs=imgsf;
for i=1:nImgs
    imgs(:,:,i)=msk.*(imgsf(:,:,i)*a+sigmaN*randn(n,n));
end;
else
   

    aS(1)=std(imgs(:))^2;
    ImgScale=sum(imgs,3)/nImgs;
    maxImgScale=max(ImgScale(:));    
    aS(2)=1.5*maxImgScale/maxRefScale; %estimate the a 


    %ctfs=ones(size(ctfs));
end;
nPkIndice=meshgrid(1:nRefs,1:nImgs)'; % to evaluate refs for each image
nImgs
%% Process images
ML.V=zeros(n,n,n,niter);
fscs=FSCorr(vol,mapx);
%if iiter==1
    alphas0=-peakAlpha:dAlpha:peakAlpha;
% else
%     alphas=
% end;
nAlphas=numel(alphas0);
IndexA=true(nAlphas,nImgs);
a=aS(2);
aI=a*ones(nImgs,1);
LlI=zeros(nImgs,1);
if b0r_reestimate==0
    b0r=0;
else
    b0r=-20/pixA/0.8; % bobbling offset 30 Angstrom
end;
%%
for iiter=1:niter
    iiter
    vesRs=vesRs-b0r;
    refs=rsMakeTemplates(angs,vol);
    
nRefs
sigmaN=sqrt(aS(1));

sigmaC=sqrt(aS(3));
b0=aS(4);


% figure(1); SetGrayscale;
% subplot(221);

%msk=fuzzymask(n,2,n*0.4,n*.1);


%r=Radius(n);
% logPC=-r.^2/(2*sigmaC^2)-log(sigmaC);
% logProbWC=repmat(logPC,[1,1,nAlphas]);

totAccumImg=zeros(n,n);
allAlphas=zeros(nAlphas,nImgs);

accumNorms=zeros(n^2,nRefs);


classSums=zeros(n^2,nRefs);
classSums1=classSums;
classSums2=classSums;
classNorms=zeros(n^2,nRefs);
classNorms1=classNorms;
classNorms2=classNorms;
accumSums=zeros(npars,1);
accumsI=zeros(npars,nImgs);
accumL=ones(nImgs,1)/nImgs;
PAlphas=zeros(numel(alphas0),nImgs);
%accumPAlphasRI=zeros(nAlphas,nRefs,nImgs);
%accumPalphasBI=zeros(nAlphas,nRefBetas,nImgs); 
tic
figure(1);SetGrayscale;subplot(111);
% if Fake
%     ImgUsed=1:nImgs;
%      writerObj = VideoWriter([num2str(iiter),'outI.avi']); % Name it.
%      writerObj.FrameRate = 10; % How many frames per second.
%      open(writerObj); 
% else
%     ImgUsed=find(si.miIndex)';
%     writerObj = VideoWriter([num2str(iiter),'toutExpI30.avi']); % Name it.
%     writerObj.FrameRate = 10; % How many frames per second.
%     open(writerObj); 
% end;
% for iImg=1:70%ImgUsed 
parfor iImg=1:nImgs
    if (ctfIndex(iImg))
    accumImages=zeros(n^2,nRefs);
    accumPars=zeros(npars,nRefs);
    logPs=-inf(nRefs,1);
    alphas=alphas0(IndexA(:,iImg));
    nAlphas=numel(alphas);
    accumPAlphas=zeros(nAlphas,1);
    %PAlphasB=zeros(nAlphas, nRefBetas,nRefGammas);
    PAlphasR=zeros(nAlphas,nRefs);
    pTransR=zeros(n*n,nAlphas,nRefs);
    funrot=fft2(rsRotateImage(imgs(:,:,iImg),alphas));

%         if ctfIndex(iImg) 
%             ctf=ctfs(:,:,ctfIndex(iImg));
%         else
%             %ctf=ctfs(:,:,1);
%             ctf=ones(n,n)*1E-9;  % zeros cannot be used
%         end;
    %ctf=ones(n,n);
    ctf=ctfs(:,:,ctfIndex(iImg));
        vesR=vesRs(iImg);
    
    nPkIndex=nPkIndice(:,iImg);
%     Loop over references
    for k=1:nRefs
        if find(nPkIndex==k)%k<=nPkIndex(k)
        beta=angs(k,2);
        logProbWC=ComputeLogProbWC2(n,sigmaC,sigmaG,vesR,y0s(iImg),alphas,dAlpha,beta,dBeta);
        %logProbWC=ComputeLogProbWC2(n,sigmaC,sigmaG,vesR+aS(4),y0s(iImg),alphas,dAlpha,beta,dBeta);
%         logProbWC=zeros(n,n,numel(alphas));
%        cref=real(ifftn(fftn(refs(:,:,k)).*ctfs(:,:,iImg)));
        cref=real(ifftn(fftn(refs(:,:,k)).*ctf)).*msk;
        [logP, accumImg, accumPar, pAlphas, pTrans]=rsOneImgOneRef5(funrot,logProbWC,aI(iImg),sigmaN,cref); 
        %loP: one element; pALphas:nAlphas pTrans: n*ns
        ctAccumImg=real(ifftn(fftn(accumImg).*ctf));
        accumImages(:,k)=ctAccumImg(:);
        logPs(k)=logP; 
        accumPar(4)=accumPar(4)/sind(beta);
        accumPars(:,k)=accumPar;
        if sum(isnan(accumPar))
            disp([num2str([iImg,k,accumPar(:)]),' accumPar']);
        end;
        %accumPAlphas=accumPAlphas+pAlphas;
        PAlphasR(:,k)=pAlphas'; % p(alpha|X,ref)      nAlphas * nRefs
        %PAlphasB(:,fix(k,nRefGammas),mod(k,nRefGammas)+1)=pAlphas'; %  each ref and alpha
        pTransR(:,:,k)=reshape(pTrans,n*n,nAlphas);  % n*n nalpha nref
        end;
    end; % k loop over refs
    maxLogP=max(logPs); %max of k elements 
    sclPk=exp(logPs-maxLogP); %k elements
    norm=sum(sclPk);  % sum_ref: p(X|ref)=p(X,ref) (ref evenly distributed)
    normPk=sclPk/norm;  % normalized p(k|Img) or p(ref_k|X)=P(X,ref_k)/sum_k(P(X,ref_k)
    normPaR=PAlphasR*diag(normPk); % nAlphas*nRefs  x nRefs*nRefs= nAlphas * nRfes
    normPa=PAlphasR*normPk; %p(alpha, X,ref)=p(alpha|X,ref)*p(X|ref)  nAlphas*nRefs  x nRefs*1  =nAlphas*1
    %(4) P(X,alpha)=sum_ref(P(X,alpha,ref))     nAlphas
    %norma=sum(normPa,2);
    %accumPAlphasRI(:,:,iImg)=normPa./repmat(norma, [1,nRefs]); 
    temp=PAlphas(:,iImg);
    temp(IndexA(:,iImg))=normPa/(sum(normPa)); % nalpha 
    PAlphas(:,iImg)=temp; % nAlphas,* nImg
    %PAlphas(IndexA(:,iImg),iImg)=norma/sum(norma); % nAlphas,* nImg
    %P(alpha|X)=P(X,alpha)/sum_alpha(P(X,alpha))	
    
    %accumPalphasBI(:,:,iImg)=sum(normPa())
%     for iBeta=1:nRefBetas
%         accumBeta(iBeta,iGamma)=sum(normPk((0+iBeta):nRefGammas:(nRef-nRefGammas+iBeta)));
%     end;

        % to get p(t|X):
       tempt=normPaR/sum(normPa(:));% nalpha*nref (6) P(alpha,ref|X)=P(X,alpha,ref)/sum_alpha,ref(P(X, alpha,ref))
        trans=reshape(pTransR,n*n,nAlphas*nRefs)*tempt(:);% (7) P(t,alpha,ref|X)=P(t|X,alpha,ref)*P(alpha,ref|X)
       % (8) p(t|X)=sum_alpha,ref(P(t,alpha,ref|x))
%        if iImg<71
%            imacs(reshape(trans,n,n));
%             colorbar;
%             title(num2str(iImg));
%             drawnow;
%        
%             frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%             writeVideo(writerObj, frame);
%           
%        end;
    nPkIndex=find(normPk>0.0001); % index of references for each image in this iteration

    if iiter==1
        nPkIndex0=zeros(nRefs,1);
        nPkIndex0(1:numel(nPkIndex))=nPkIndex;
        %nPkIndice(:,iImg)=0;
        nPkIndice(:,iImg)=nPkIndex0;
          
        %if (max(PAlphas)./(min(PAlphas)+1E-10))>0.001
    end;
    aI(iImg)=accumPars(2,:)*normPk;
    %accumL=norm*exp(maxLogP)+accumL; 
    LlI(iImg)=log(norm)+maxLogP; % per image
    accumSums=accumSums+accumPars*normPk;
    accumsI(:,iImg)=accumPars*normPk;
    pixNormPk=repmat(normPk',n^2,1);
    
    classSums1=classSums1+accumImages.*repmat(normPk',n^2,1)*mod(iImg,2);  % npixels x nrefs
    classNorms1=classNorms1+ctf(:).^2*normPk'*mod(iImg,2);
    classSums2=classSums2+accumImages.*repmat(normPk',n^2,1)*(~mod(iImg,2));
    classNorms2=classNorms2+ctf(:).^2*normPk'*(~mod(iImg,2));
% %     Show the best aligned image
%     subplot(221);
%     [mxVal mxK]=max(normPk);  % show the reference with the best alpha
%     imacs(reshape(accumImages(:,mxK),n,n));
%     title(num2str([iImg mxK]));
%     
%     subplot(222);
%     plot(alphas,accumPAlphas/iImg);
%     ylabel('pAlphas');
%     xlabel('Alpha');
%     drawnow;
    end;  % if ctfIndex
%     if iImg==70
%         close(writerObj); % Saves the movie.  
%     end;
end; % iImg loop over images
toc
%%  Display the results

    if iiter==1
        IndexA=PAlphas>0.001;
    end;
    maxlogPI=max(LlI);
    accumL=sum(exp(LlI-maxlogPI)); % all imges
accumLl=log(accumL)+maxlogPI
figure(1);
%disp(accumSums/nImgs);  % accumulators
subplot(222);
plot(accumsI(2,:)); title('Amplitudes');
xlabel('Image no.');
set(gca,'fontsize',14);

figure(2);
SetGrayscale;
clf;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%need editting for
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%fsc real
classSums=classSums1+classSums2;
classNorms=classNorms1+classNorms2;
clsSum=reshape(classSums,n,n,nRefs);
clsNorm=reshape(classNorms,n,n,nRefs);
k=.2;
classMeans=real(ifft2(fft2(clsSum)./(k+clsNorm)));
ixs=.5:nRefGammas+.5;
iys=.5:nRefBetas+.5;
imacs(ixs,iys,ImageArray(classMeans,0,n,1,nRefGammas,nRefBetas));
set(gca,'xtick',1:nRefGammas,'XTickLabel',num2str(refGammas'));
xlabel('Gamma');
set(gca,'ytick',1:nRefBetas,'YTickLabel',num2str(refBetas'));
ylabel('Beta');

figure(1);
subplot(223);
semilogy(sum(classNorms),'-');
axis([0 nRefs 1e-8 inf]);
set(gca,'fontsize',14);
% hold on; % Mark the references of the original images
% semilogy(refPtrs,classNorms(refPtrs),'r+');
% hold off;
xlabel('Reference no.');
ylabel('Class norm');

subplot(224);
[mxVal mxK]=max(sum(classNorms));
bestImg=reshape(classSums(:,mxK),n,n);
imacs(bestImg);
title('Strongest class mean');
set(gca,'fontsize',14);
aS=accumSums/sum(logical(find(ctfIndex>0)));
%aS=accumSums/nImgs;
aS(1)=aS(1)*n^2/(msk(:)'*msk(:));
% aS(3)=aS(3)-aS(4)^2;
%theoret=[sigmaN^2 a sigmaC^2 b0 0];
if ~Fake
    theor_a=aS(2)*sVesicle/mean(sVesicle);
    ctfIndex(logical((aI./theor_a>1.8)+(aI./theor_a<0.55)))=0;
else
    ctfIndex(logical((aI<1/10*aS(2))+(aI>10*aS(2))))=0; % neglect imgs with absurd a
end;
textCells={'noise var  '
           'a          '
           'click var  '
           'b0         '
           '           '};
for i=1:numel(accumSums)
    disp([textCells{i} num2str([aS(i) theoret(i)])]);
end;
%%
% clsNorm=repmat(ifftshift(fuzzymask(n,2,0.4*n,.05*n)),[1 1 nRefs]);
% rclsNorm=zeros(n,n,nRefs);
% clsSum=refs;
symmetry=2;
k=1;% for k2camera 0.1 is too small(noise start to incease) 0.4-1 looks the same
figure(2+iiter);
if ~(GoldFSC||iiter<niter)
    % FSC
    %fsc=FSCbyhalves(classSums1,classSums2,classNorms1,classNorms2,n,nRefs,angs);

    %[fsc,mapUM]=VFSCbyhalvesM(k,symmetry,classSums1,classSums2,classNorms1,classNorms2,n,nRefs,angs,refMsks);
    [fsc,mapUM]=VFSCbyhalves(k,symmetry,classSums1,classSums2,classNorms1,classNorms2,n,nRefs,angs);
else
    rclsNorm=zeros(n,n,nRefs);
    for i=1:nRefs
        rclsNorm(:,:,i)=fftshift(ifft2(clsNorm(:,:,i)));
    end;
    [reconVol normVol]=rsDoReconstructionSymmetry(clsSum,rclsNorm,angs,symmetry);
    % reconstruction
    mapUM=rsNormalizeReconstruction(reconVol,normVol,k);
    fsc=FSCorr(vol,mapx);  % Compare with crystal structure   
end;

mapUM_f=SharpFilt(mapUM,0.5,0.025);
map=mapUM_f.*Volmsk;
vol=map*std(vol(:))/std(map(:));  % scale
ShowSections(vol);
fscs=[fscs,fsc];
subplot(339);
plot(fscs);
set(gca,'xtick',1:3:n/2,'XTickLabel',round(n*pixA./(1:3:n/2)'));
xlabel('resolution (Å)');
line(pixA*n./[30,30],[0,1]);
line(pixA*n./[25,25],[0,1]);
line(pixA*n./[20,20],[0,1]);
line([0,n/2],[0.5,0.5]);
line([0,n/2],[0.143,0.143]);
%% store the values
ML.V(:,:,:,iiter)=vol;
ML.sigmaN2(iiter)=aS(1);
ML.a(iiter)=aS(2);
ML.sigmaC2(iiter)=aS(3);
ML.b0(iiter)=aS(4);
ML.L(iiter)=accumLl;
ML.fscs=fscs;
ML.pixA=pixA;
if b0r_reestimate %% && abs(aS(4))>3
    b0r=aS(4);
else
    b0r=0;
end;
save(['Results' basename '.mat'],'ML');
end;
disp('Done.');

vol(:,:,membrane-round(40/pixA):membrane)=vol(:,:,membrane-round(40/pixA):membrane)*3/2;
outname=[basename 'm.mrc'];
WriteMRC(vol,pixA,outname);%'300EFm_b10r15.mrc');
WriteMRC(mapUM,pixA,[basename '_UM.mrc']);
%%
% if GoldFSC
%     FSCgold_user;
% end


