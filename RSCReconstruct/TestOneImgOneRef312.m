%function TestOneImgOneRef3
% version 3 includes ctfs
% Test the likelihood calculation of 2D images
%version 9 correct ctfs format, and make a as a vectors to filter out
%absurds elements
% V10 change ctfs=... (abs(CTF(...))) for fake images
% version 11 is for new datasets
% VERSION 17 EXCLUDE BAD images by low pass filter
% V18 can makes normal fscs  line 440 and 441 to switch and make symmetry
% =2.  Now put in dataset with second exposure, this stack has been 
%eliminated.so  line 165 %ctfIndex(exc)=0;
% Jun 26 mask the volume, change fscbyhalves to vFSCbyhalves, and use the
% right reconstruction in symmetry.

% V28 apply msk before reconstruction good Jun30, corrected the mask for 2D July2013 
% V30 is modified form V28, change (1) P(X,ref,alphas) line 311 (2)normPa
% Line 312  (3) accumL scale to maxlog
% V31 introduces nmodel=2 (ref to V23)

Fake=1;


if Fake
    n=40;
    disp('This is a simulation');
else
    n=48;
    disp('Operated on experimental images')
end;
ds=2;
dAlpha=5;
peakAlpha=10;
dBeta=5;
a=.01;
vesR=100/ds;
nrep=6;

niter=3;

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
nmodel=2;
aS(2)=0.4;  %initial guess of parameters
aS(1)=0.8;
%sigmaR=5;
aS(3)=8;
%sigmaG=2;
aS(4)=2;
%aS(4)=0.5;
aS=theoret;

load /Users/fred/aEMCodeRepository/RSC/3KG2RotMap2.9A.mat
%load /home/yc323/RSC-ML2/3KG2map5.8A.mat
%load /Users/yi/Installers/EMData/RSC201304/RSC/3KG2map5.8A.mat
% load ../RSC2/3KG2RotMap2.9A.mat
% load ../RSC/3KG2RotMap2.9A.mat
pixA0=2.9;  % model pixel size
nm=size(map,1);ds=6.97/(pixA0); %6.97 3.486is the pixel size for the stack
mapx=Crop(DownsampleGeneral(map,round(nm/ds)),n);
pixA=pixA0*ds;

membrane=round(n/2-35/pixA);
%mapx=Crop(Downsample(map,nm/ds),n);
mapx(:,:,1:membrane)=mapx(:,:,1:membrane)*2/3;
%mapx(:,:,1+4:n)=mapx(:,:,1:n-4);
%mapx(:,:,1:n-3)=mapx(:,:,1+3:n);
mapx1=mapx;
for i=membrane:n
    mapx1(:,:,i)=DownsampleGeneral(mapx(:,:,i),n,1.2);
end;
mapx2=mapx1;
for i=1:n
    tempmap=DownsampleGeneral(squeeze(mapx1(i,:,:)),n,1/1.2);
    mapx2(i,:,membrane:n)=tempmap(:,membrane:n);
end;

RefScale=sum(mapx,3);
maxRefScale=max(RefScale(:));
fc=1/(60/pixA);

vols=zeros(n,n,n,nmodel);
vols(:,:,:,1)=SharpFilt(mapx,fc,0);
vols(:,:,:,2)=SharpFilt(mapx2,fc,0);

Volmsk=(GaussFilt(mapx,fc/2)>0.01*std(mapx(:)));
vol=SharpFilt(mapx,fc,0);

figure(2);
SetGrayscale;

%% Make references
angs=[];
k=0;    
%refBetas=10:10:170;
refBetas=30:10:150;
refGammas=0:15:165;
%refBetas=20:20:160;
%refGammas=0:30:165;
nRefBetas=numel(refBetas);
nRefGammas=numel(refGammas);
for beta=refBetas
    for gamma=refGammas;  % There are 6 gammas
        k=k+1;
        angs(k,:)=[0 beta gamma];
    end;
end;
%refMsks=rsMakeTemplates(angs,Volmsk);  this is not bia-free
refs=GaussFilt(rsMakeTemplates(angs,mapx),0.05,1);
refMsks0=refs>0.01*std(refs(:));
refMsks=GaussFilt(refMsks0,0.05,1);
nRefs0=k;
nRefs=k*nmodel;
allBetas=angs(:,2);
betaRange=[min(allBetas,0) allBetas+5];
priorPKs=mean(sind(betaRange),2).*(betaRange(:,2)-betaRange(:,1))*pi^2/180^2;
%nRefBetas, % add  along each beta, to evaulate alphas fro all the images 

%refs=rsMakeTemplates(angs,mapx);

% figure(3);
% ImagicDisplay1(refs,2);

%% Create the fake images or load images
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
    %load sq03_10027tsi.mat;% 3124 particles with merged ctf
    %imgs0=ReadMRC('sq03_10027stall.mrc');
    load ../RSC-ML2/sq03_10027w010e1tsi.mat;% 3124 particles with second exposure
    imgs0=ReadMRC('../RSC-ML2/sq03_10027w010e1stall.mrc');
    
    y0s=si.yClick;
    nImgs=numel(y0s); 
    n0=size(imgs0,1); 
    if n==n0
        imgs=imgs0;
    end;
%    pixA0=si.pixA;
%     ntemp=round(n0*pixA0/pixA); % for downsampling
%     tempimg=zeros(ntemp,ntemp,nImgs);
%     for i=1:nImgs
%         tempimg(:,:,i)=DownsampleGeneral(imgs0(:,:,i),ntemp);
%     end;
%     imgs=Crop(tempimg,n,1);
    %imgs=imgs(n0/2-n/2+1:n0/2+n/2,n0/2-n/2+1:n0/2+n/2,:);
    
    ctfs=ifftshift(si.ctfs);
    vesRs=si.rVesicle;
    ctfIndex=si.miIndex;
%% exclude particles
    d0=160;  % angstroms min freqency
d1=80; % angstroms
    sps=RadialPowerSpectrum(imgs);
        fmin=floor(si.pixA*n/d0);
        fmax=ceil(si.pixA*n/d1);
        lfv=(fmin:fmax)*sps(fmin:fmax,:);
        sigmap=std(lfv);
        threshSD=2.5;
    exc=lfv>threshSD*sigmap;
    %ctfIndex(exc)=0;
    aS(1)=std(imgs(:))^2;
    ImgScale=sum(imgs,3)/nImgs;
    maxImgScale=max(ImgScale(:));    
    aS(2)=1.5*maxImgScale/maxRefScale; %estimate the a 
    msk=fuzzymask(n,2,0.40*n,.1*n);
    for i=1:nImgs
        imgs(:,:,i)=msk.*(imgs(:,:,i));
    end;

    %ctfs=ones(size(ctfs));
end;
nPkIndice=meshgrid(1:nRefs,1:nImgs)'; % to evaluate refs for each image
nImgs
%% Process images
%ML.V=zeros(n,n,n,niter);
%fscs=FSCorr(vol,mapx);
ML.V{nmodel}=zeros(n,n,n,niter);
fscs=FSCorr(vols(:,:,:,1),mapx);
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

fmodel=ones(nRefs,1)/nmodel;
%%
for iiter=1:niter
    iiter
    LlI1=zeros(nImgs,1);
    LlI2=zeros(nImgs,1);% for different models
    for imodel=1:nmodel
        vol=vols(:,:,:,imodel);
        refs0=rsMakeTemplates(angs,vol);
        refs(:,:,nRefs0*(imodel-1)+1:nRefs0*imodel)=refs0;
    end;
    
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
 %for iImg=1:70%ImgUsed 
% parfor iImg=1:nImgs
for iImg=1:nImgs
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
        if mod(k,nRefs0)
                            beta=angs(mod(k,nRefs0),2);
                        else
                            beta=angs(nRefs0,2);
                        end
        logProbWC=ComputeLogProbWC2(n,sigmaC,sigmaG,vesR,y0s(iImg),alphas,dAlpha,beta,dBeta);
        %logProbWC=ComputeLogProbWC2(n,sigmaC,sigmaG,vesR+aS(4),y0s(iImg),alphas,dAlpha,beta,dBeta);
%         logProbWC=zeros(n,n,numel(alphas));
%        cref=real(ifftn(fftn(refs(:,:,k)).*ctfs(:,:,iImg)));
        cref=real(ifftn(fftn(refs(:,:,k)).*ctf)).*msk;
        [logP accumImg accumPar pAlphas pTrans]=rsOneImgOneRef5(funrot,logProbWC,aI(iImg),sigmaN,cref); 
        %loP: one element; pALphas:nAlphas pTrans: n*ns
        ctAccumImg=real(ifftn(fftn(accumImg).*ctf));
        accumImages(:,k)=ctAccumImg(:);
        logPs(k)=logP; 
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
    logPs=logPs+log(fmodel);
    maxLogP=max(logPs); %max of k elements 
    sclPk=exp(logPs-maxLogP); %k elements
% %     [norm1,norm2]=sum(reshape(sclPk,nRefs0,nmodel),2);
    norm1=sum(sclPk(1:nRefs0));
    norm2=sum(sclPk(1+nRefs0:nRefs));
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
    LlI1(iImg)=log(norm1)+maxLogP;
    LlI2(iImg)=log(norm2)+maxLogP;
    
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
    accumL1=sum(exp(LlI1-maxlogPI));
    accumL2=sum(exp(LlI2-maxlogPI));
    fmodelI=[accumL1,accumL2]/(accumL1+accumL2)
    fmodel=reshape(repmat(fmodelI,nRefs0,nmodel),nRefs,1);
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
imacs(ixs,iys,ImageArray(classMeans,0,n,1,nRefGammas,nmodel*nRefBetas));
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
    theor_a=aS(2)*si.sVesicle/mean(si.sVesicle);
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
rclsNorm=zeros(n,n,nRefs);
for i=1:nRefs
    rclsNorm(:,:,i)=fftshift(ifft2(clsNorm(:,:,i)));
end;
symmetry=2;
%%
%[reconVol normVol]=rsDoReconstructionSymmetry(clsSum,rclsNorm,angs,symmetry);
for imodel=1:nmodel
figure(2+(iiter-1)*nmodel+imodel);
% reconstruction
k=1;
%mapUM=rsNormalizeReconstruction(reconVol,normVol,k);


% FSC
%fsc=FSCbyhalves(classSums1,classSums2,classNorms1,classNorms2,n,nRefs,angs);
%fsc=FSCorr(vol,mapx);  % Compare with crystal structure

    range=nRefs0*(imodel-1)+1:nRefs0*imodel;
[fsc,mapUM]=VFSCbyhalvesM(k,symmetry,classSums1(:,range),classSums2(:,range),classNorms1(:,range),classNorms2(:,range),n,nRefs0,angs,refMsks);
map=mapUM.*Volmsk;
vol=map*std(vol(:))/std(map(:));  % scale
ShowSections(vol);
fscs=[fscs,fsc];
subplot(339);
%plot(fscs);
plot(fscs(:,1+imodel:nmodel:nfsc));
set(gca,'xtick',1:3:n/2,'XTickLabel',round(n*pixA./(1:3:n/2)'));
xlabel('resolution (Å)');
line(pixA*n./[30,30],[0,1]);
line(pixA*n./[25,25],[0,1]);
line(pixA*n./[20,20],[0,1]);
line([0,n/2],[0.5,0.5]);
ML.V{imodel}(:,:,:,iiter)=vol;
end;
%% store the values

ML.sigmaN2(iiter)=aS(1);
ML.a(iiter)=aS(2);
ML.sigmaC2(iiter)=aS(3);
ML.b0(iiter)=aS(4);
Ml.L(iiter)=accumLl;
end;
disp('Done.');
save('ResultsOldData10Inter.mat','ML');


