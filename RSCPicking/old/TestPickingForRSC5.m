% TestPickingForRSC5(m,mi,map,mapPixA)
% Generic inverse filter params
f0=.02;  % 4th order hp
% f0=.05;
ex=4;
f1=.0015; % Lorentizian
% f1=.004; % Lorentizian
a1=.4;

ds1=2;  % extra downsampling.
simulateImage=0;


if simulateImage
    TestFakeMicrograph;
    minS=.055;
    maxS=.075;
    origImg=simImage;
else
    
    basePath='/Users/fred/matlabWork/Yunhui2/VesicleFinderData/120711/';
    miName='Info/004_sq02_1_04mi.mat';
    %     miName='InfoTestRefinement/007_sq02_1_07smi.mat';
    %     baseName='004_sq02_1_04';
    load([basePath miName]);
    amps=mi.vesicle.s;
    %     minS=min(amps)
    %     maxS=max(amps)
    minS=.055;
    maxS=.075;
    mi.vesicle.ok=(amps>minS) & (amps<maxS);
    
    baseName=mi.baseFilename;
    imgName=['Merged/' baseName 'm.mrc'];
    origImg=ReadEMFile([basePath imgName]);
    
    % basePath='/Volumes/TetraData/EMWork/Hideki/120811/AMPA_R_lipo_washed/';
    % baseName='sq02_1_01';
    % miName=['Info/' baseName 'mi.mat'];
    % imgName=['Merged/' baseName 'm.mrc'];
    
    % minS=.085;
    % maxS=.105;
    % load([basePath miName]);
    
end;

n=size(origImg);

if ds1>1  % further downsampling
    n=n/ds1;
    m0=Downsample(origImg,n);
else
    m0=origImg;
end;

ds=mi.imageSize(1)/n(1);

if simulateImage
    m1=Downsample(mParticles,n);
else
    disp('Making model vesicles');
    ves=meMakeModelVesicles(mi,n);
    m1=m0-ves;
end;

disp('High-pass filtering');
% perform high-pass (inverse) filtering
f=RadiusNorm(n)/(mi.pixA*ds);  % frequencies for evaluating the CTF
hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
hinv(f==0)=0;
m=real(ifftn(fftn(m1).*ifftshift(hinv)));

figure(1); clf;
SetGrayscale;
imac(imscale(GaussFilt(m,.2),256,.001));
drawnow;

%% Create the templates
load /Volumes/TetraData/Structures/AMPAR/3KG2map5.8AGoodScale.mat
nt0=size(map,1)/ds1;
if ds1>1
    map=Downsample(map,nt0);
end;
nt=ceil(nt0/8)*8;  % has to be a multiple of 8 for gridding fcns.
map=Crop(map,nt);
nt

map=map*5.8*ds1;  % approx amplitude correction (V-A scaling)
ds=2*ds1;
membraneOffset=-24/ds;  % downsampled map by 2.
nAlpha=32; % about 10 degrees
nBeta=12;  % even is best.  Number on hemisphere.
nGamma=8;
symmetry=2;
gammaStep=360/(symmetry*nGamma);
nterms=26;
nterms=12;

% nterms=200;

[hemiAngles angleInds]=rsListHemisphereAngles(nAlpha, nBeta);
nHemiAngles=size(hemiAngles,1);
nHemi=2;  % both hemispheres

angleList=zeros(nHemi,nGamma,nHemiAngles,3);
for j=1:nGamma;
    gamma=(j-1)*gammaStep;
    for k=1:nHemiAngles
        angleList(1,j,k,:)=[hemiAngles(k,:) gamma];
        angleList(2,j,k,:)=[[0 180]-hemiAngles(k,:) gamma];
    end;
end;
nAngles=numel(angleList)/3;

% angle list is of size
% (nHemi x nGamma x nHemiAngles, 3)  where nHemi=2 is the number of hemispheres.
% Note that the beta angles alternate such that, if betastep=1, they are
% are (0, 180) for each gamma, then (1, 179) for each gamma, up to 89, 91.
% that is, the same projected position is described twice.

%%

disp(['Making ' num2str(nAngles) ' templates']);

tic
% allTemplates=rsMakeTemplatesQuick(angleList,map);
allTemplates=rsMakeTemplates(reshape(angleList,nAngles,3),map);
toc

% angleList=reshape(angleList,nHemiAngles,2,nGamma,3);
allTemplates=reshape(allTemplates,nt,nt,nHemi,nGamma,nHemiAngles);


%% Filter the templates according to the CTF
[nt nt nHemi nGamma nHemiAngles]=size(allTemplates);
nAngles=nHemi*nGamma*nHemiAngles;
ne=NextNiceNumber(nt*1.2);  % increase the size to allow CTF rings
ctf=meGetEffectiveCTF(mi,ne,ds);  % put in dqe, pw filter.
% evaluate a generic inverse filter
f=RadiusNorm(ne)/(mi.pixA*ds);  % frequencies for evaluating the CTF
hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
hinv(f==0)=0;
% Pad the templates to avoid ctf artifacts
% eigenSet=rsMakeEigenrefs(allTemplates,nterms,ctf.*hinv);
xTemplates=Crop(reshape(allTemplates,nt,nt,nAngles),ne,1);
nim=size(xTemplates,3);
H=ifftshift(hinv.*ctf);
msk=fuzzymask(ne,2,0.45*ne,.1*ne);
for i=1:nim
    xTemplates(:,:,i)=real(ifftn(fftn(xTemplates(:,:,i)).*H)).*msk;
end;
xTemplates=reshape(xTemplates,ne,ne,nHemi,nGamma,nHemiAngles);
ntstr=num2str(nt);
nestr=num2str(ne);
disp(['Templates expanded from ' ntstr 'x' ntstr ' to ' nestr 'x' nestr]);

%% Make the eigenreferences
eigenSet=rsMakeEigenrefs(xTemplates,nterms);
figure(2);

plot(1-eigenSet.termVar);
% figure(1);
% ImagicDisplay(eigenSet.imgs,2);


%%  % Show the templates and the reconstructions
timgs=reshape(eigenSet.imgs,ne*ne,nterms);
nG=1;
nH=2;
nAngs=nHemiAngles;
rImg=single(zeros(ne,ne,2*nH,nG,nAngs));
for k=1:nHemiAngles
    for j=1:nG
        for i=1:nH
            %             rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne).*squeeze(eigenSet.ampList(i,j,k));
            rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne);
            rImg(:,:,2*i,j,k)=xTemplates(:,:,i,j,k);
        end;
    end;
end;
% figure(1);
% ImagicDisplay2(rImg,2);


%%  % Make a complete set of reconstructions
%   and make the average power spectrum
nTotal=nAngs*nGamma*nHemi;
timgs=reshape(eigenSet.imgs,ne*ne,nterms);
rImg=single(zeros(ne,ne,nTotal));
vL=reshape(eigenSet.vList,nterms,nTotal);
sp=zeros(ne,ne);
sp0=zeros(ne,ne);
for i=1:nTotal
    %     img=reshape(timgs*vL(:,i),ne,ne)*eigenSet.ampList(i);
    img=reshape(timgs*vL(:,i),ne,ne);
    sp0=sp0+abs(fftn(xTemplates(:,:,i))).^2;
    sp=sp+abs(fftn(img)).^2;
    rImg(:,:,i)=img;
end;
spr=Radial(fftshift(sp/nTotal));
spr0=Radial(fftshift(sp0/nTotal));
figure(2);
semilogy([spr0 spr]);
figure(3);
plot([cumsum(spr0) cumsum(spr)]);
% ImagicDisplay2(rImg);



%%  Evaluate the cc for all vesicles
figure(1); clf; SetGrayscale;
imac(imscale(GaussHP(GaussFilt(origImg,.1),.005),256,.003));
namp=1.5;
% namp=0.5;
if simulateImage
    m2=m+namp*randn(n);
else
    m2=m;
end;

mic=mi;  % copy
% mic.vesicle.ok=(mic.vesicle.s>minS) & (mic.vesicle.s<maxS);
% % mic.vesicle.s=mic.vesicle.s*0+.05;  % make the amplitudes equal.
vesInds=find(mic.vesicle.ok);

disp(['Total vesicles: ' num2str(numel(vesInds))]);
mxCC=zeros(n);
mxVar=zeros(n);
mxTempInds=zeros(n);
mxVesInds=zeros(n);
figure(3);
SetGrayscale;
partRadius=nt/3;


for i=vesInds
    % for i=1:16
    %     Compute the eigen-CC for one vesicle
    [mxVals mxInds var cImageX iptrs]=rsVesicleCorrelation(-m2,mic,i,membraneOffset,partRadius,angleInds,eigenSet);
    %     Merge it into the mxCC and mxGamma images
    q=mxVals>mxCC;
    mxCC(q)=mxVals(q);
    mxVesInds(q)=i;  % the related vesicle
    mxTempInds(q)=mxInds(q);
    q2=var>mxVar;
    mxVar(q2)=var(q2);
    if mod(i,10)==0
        imacs(mxCC);
        title(i);
        drawnow;
    end;
end;
mxCC8=uint8(mxCC*100);

%
if simulateImage
    nparts=8;
    % For the first row of vesicles/particles
    iTemplates=zeros(nparts,1);
    iCCVal=zeros(nparts,1);
    imxs=zeros(nparts,1);
    imys=zeros(nparts,1);
    mCopy=m;
    nl=48;  % size of local image
    mskl=fuzzymask(nl,2,17,2);
    % eigenRefs=Crop(eigenSet.imgs,nl,1);  % pad the images
    % eR=reshape(eigenRefs,nl^2,nterms);
    
    for ix=1:nparts  % first row
        %     imx=imx0+(ix-1)*64;
        ix1=(ix-1)*64;
        [v imx imy]=max2d(mxCC(ix1+1:ix1+64,1:350));  % Get the vicinity of the particle.
        imx=imx+ix1;
        
        imxs(ix)=imx;
        imys(ix)=imy;
        itempl=mxTempInds(imx,imy);
        iTemplates(ix)=itempl;
        itempind=max(1,itempl);
        iCCVal(ix)=mxCC(imx,imy);
        theTemplate=xTemplates(:,:,itempind);
        temp=ExtractImage(theTemplate,[imx imy+50],n,1);
        mCopy=mCopy-iCCVal(ix)*mi.vesicle.s(ix)*temp;  % add it to the image
        
        %    Compute the residual variance
        % First, create a fake vesicle image
        resRef=xTemplates(:,:,itempind)*mic.vesicle.s(ix);
        mx=ExtractImage(resRef,[imx imy],n,1);
        r=mi.vesicle.r(ix)/ds-membraneOffset;
        ctr=[mi.vesicle.x(ix) mi.vesicle.y(ix)]/ds+1;  % shift to 1-based coords
        [txVals txInds tvar tImageX tptrs]=...
            rsVesicleCorrelation(mx,mic,ix,membraneOffset,partRadius,angleInds,eigenSet);
        localCC=ExtractImage(mxCC,[imx imy],nl).*mskl;
        localModel=ExtractImage(txVals,[imx imy],nl).*mskl;
        residCC=(localCC-localModel);
        residVar=residCC(:)'*residCC(:);
        ccVar=localCC(:)'*localCC(:);
        [ix ccVar residVar]
        
    end
    %
    iAngles=floor((iTemplates-1)/16)+1;
    % iAngles=min(iAngles,size(hemiAngles,1));
    gma=ceil(mod(iTemplates-1,16)/2)*gammaStep;
    aList=reshape(angleList,nAngles,3);
    [imxs imys iCCVal iTemplates mod(iTemplates-1,16)+1 aList(iTemplates,:)]
    imacs(mCopy)
end;
%%
% if simulateImage
%     localCC=mxCC(imx-16:imx+15,imy-16:imy+15);
%     localTm=mxTempInds(imx-16:imx+15,imy-16:imy+15)/32;
%     imacs(localCC);
%     ccMaximum=mxCC(imx,imy)
%     ccTemplate=mxTempInds(imx,imy+1);
%     ccBestAng=floor(ccTemplate/32)
%     ccBestRef=mod(ccTemplate,32)
%     
%     allAngles=reshape(angleList,nAngles,3);
%     
%     
%     %%
%     % localCCs=cImage(:,ilx-16:ilx+15,ily-16:ily+15);
%     % plot(localCCs(:,:,17)');
% end;
%%
q=squeeze(max(localCCs));
imacs(q);
% [val ind]=max(cImage(:,ilx,ily))
% thePointer=iptrs(ilx,ily,1)



%%
% save([basePath 'Info/' baseName 'ptmp.mat'],'mxCC8');
% imwrite(rot90(mxCC8),[basePath 'Info/' baseName 'pcc.tif']);
imwrite(rot90(mxCC8),[basePath baseName 'pcc.tif']);

% %% Make variance map
% rParticle=10;  % pixels
% disc=ifftshift(fuzzymask(n,2,rParticle,rParticle/10));
% disc=disc/sum(disc(:));  % unit overall amplitude
% filtVar=real(ifftn(fftn(mxVar).*fftn(disc)));
% filtCC=real(ifftn(fftn(mxCC).*fftn(disc)));

%% Search for CC peaks

blankRadius=16;
minAmp=.6;
maxAmp=1.5;
mxCC2=mxCC;
nFound=0;
maxFound=100;
iter=0;
maxIters=50;

threshVar=10;
mskl=fuzzymask(nl,2,17,2);
nomask=1;
[amp ix iy]=max2d(mxCC2);
coords=[];
amps=zeros(2,1);
figure(2);
while amp>minAmp && nFound < maxFound && iter < maxIters
    iter=iter+1
    if amp <= maxAmp && amp > minAmp  % valid peak
        itempind=max(1,mxTempInds(ix,iy));
        
        %    Compute the residual variance
        % First, get the template
        iv=mxVesInds(ix,iy);  % Look up the vesicle index
        [iv ix iy]
%         make the model image assuming amp=1.
        resRef=xTemplates(:,:,itempind)*mic.vesicle.s(iv);
        subplot(222);
        imacs(resRef);
        
        %     Create a full fake image and correlate it.
        mx=ExtractImage(resRef,[ix iy],n,1);
        r=mi.vesicle.r(iv)/ds-membraneOffset;
        ctr=[mi.vesicle.x(iv) mi.vesicle.y(iv)]/ds+1;  % shift to 1-based coords
        [txVals txInds]=rsVesicleCorrelation(mx,mic,iv,membraneOffset,...
            partRadius,angleInds,eigenSet);
        localCC=ExtractImage(mxCC,[ix iy],nl).*mskl;
        nCCPoints=sum(localCC(:)~=0);

        subplot(221); imacs(localCC);
        localModel=ExtractImage(txVals,[ix iy],nl).*mskl;
        residCC=(localCC-localModel);
        subplot(223)
        imacs(residCC);
%         imacs(localModel);
        residVar=residCC(:)'*residCC(:);
        ccVar=localCC(:)'*localCC(:);
        varStats=[amp ccVar residVar]*sum(mskl(:))/nCCPoints
        if residVar<threshVar
            nFound=nFound+1;
            coords(nFound,:)=[ix iy]; %-n/2-1;
            amps(nFound)=amp;
            found=[nFound amp ix iy]
        end;
        mxCC2=mxCC2.*(1-fuzzymask(n,2,blankRadius,2,[ix iy]));
        subplot(224);
        imacs(mxCC2);
        title(nFound);
        drawnow;
        [amp ix iy]=max2d(mxCC2);
    end;
end;
