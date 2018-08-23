% TestPickingForRSC4(m,mi,map,mapPixA)
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

disp('Making model vesicles');
ves=meMakeModelVesicles(mi,n);
m1=m0-ves;

disp('High-pass filtering');
% perform high-pass (inverse) filtering
f=RadiusNorm(n)/(mi.pixA*ds);  % frequencies for evaluating the CTF
hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
hinv(f==0)=0;
m=real(ifftn(fftn(m1).*ifftshift(hinv)));

% m=-mParticles+2*randn(n);

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
% nterms=200;

[hemiAngles angleInds]=rsListHemisphereAngles(nAlpha, nBeta);
nHemiAngles=size(hemiAngles,1);

% make two copies, one for upper and one for lower hemisphere.
sphereAngles=[hemiAngles
    repmat([0 180],nHemiAngles,1)-hemiAngles];
angleList=repmat(sphereAngles,nGamma,1);
% increment gamma each 2*nHemiAngles
gammas=repmat((0:nGamma-1)*gammaStep,2*nHemiAngles,1);
angleList(:,3)=gammas(:);  % gamma angles are most slowly varying.
nAngles=size(angleList,1);
% angle list is of size
% (nHemiAngles x nH x nGamma, 3)  where nH=2 is the number of hemispheres.
% Note that the nHemiAngles goes from beta=0 to 89, then beta =  180 to 91;
% that is, the same projected position is described twice.  Then all this
% is repeated nGamma times, with the gamma angle incremented each time.

%%
disp(['Making ' num2str(nAngles) ' templates']);

tic
% allTemplates=rsMakeTemplatesQuick(angleList,map);
allTemplates=rsMakeTemplates(angleList,map);
toc

angleList=reshape(angleList,nHemiAngles,2,nGamma,3);
allTemplates=reshape(allTemplates,nt,nt,nHemiAngles,2,nGamma);


%% Filter the templates according to the CTF
[nt nt nAngs nHemi nGamma]=size(allTemplates);
ne=NextNiceNumber(nt*1.2);  % increase the size to allow CTF rings
ctf=meGetEffectiveCTF(mi,ne,ds);  % put in dqe, pw filter.
% evaluate a generic inverse filter
f=RadiusNorm(ne)/(mi.pixA*ds);  % frequencies for evaluating the CTF
hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
hinv(f==0)=0;
% Pad the templates to avoid ctf artifacts
% eigenSet=rsMakeEigenrefs(allTemplates,nterms,ctf.*hinv);
xTemplates=Crop(reshape(allTemplates,nt,nt,nAngs*nHemi*nGamma),ne,1);
nim=size(xTemplates,3);
H=ifftshift(hinv.*ctf);
msk=fuzzymask(ne,2,0.45*ne,.1*ne);
for i=1:nim
    xTemplates(:,:,i)=real(ifftn(fftn(xTemplates(:,:,i)).*H)).*msk;
end;
xTemplates=reshape(xTemplates,ne,ne,nAngs,nHemi,nGamma);
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
nH=1;
rImg=single(zeros(ne,ne,2*nAngs,nG,nH));
for k=1:nG
    for j=1:nH
        for i=1:nAngs
            rImg(:,:,2*i-1,j,k)=reshape(timgs*eigenSet.vList(:,i,j,k),ne,ne).*squeeze(eigenSet.ampList(i,j,k));
            rImg(:,:,2*i,j,k)=xTemplates(:,:,i,j,k);
        end;
    end;
end;
figure(1);
ImagicDisplay2(rImg);


%%  % Make a complete set of reconstructions
%   and make the average power spectrum
nTotal=nAngs*nGamma*nHemi;
timgs=reshape(eigenSet.imgs,ne*ne,nterms);
rImg=single(zeros(ne,ne,nTotal));
vL=reshape(eigenSet.vList,nterms,nTotal);
sp=zeros(ne,ne);
sp0=zeros(ne,ne);
for i=1:nTotal
    img=reshape(timgs*vL(:,i),ne,ne)*eigenSet.ampList(i);
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


mic=mi;  % copy
% mic.vesicle.ok=(mic.vesicle.s>minS) & (mic.vesicle.s<maxS);
% % mic.vesicle.s=mic.vesicle.s*0+.05;  % make the amplitudes equal.
vesInds=find(mic.vesicle.ok);
disp(['Total vesicles: ' num2str(numel(vesInds))]);
mxCC=zeros(n);
mxVar=zeros(n);
mxGamma=zeros(n);
mxAngles=zeros(n);

figure(3);
SetGrayscale;
partRadius=nt/3;

for i=vesInds
    %     Compute the eigen-CC for one vesicle
    [mxVals mxInds mxAngs var]=rsVesicleCorrelation(-m,mic,i,membraneOffset,partRadius,angleInds,eigenSet);
    %     Merge it into the mxCC and mxGamma images
    q=mxVals>mxCC;
    mxCC(q)=mxVals(q);
    mxGamma(q)=mxInds(q);
    mxAngles(q)=mxAngs(q);
    q2=var>mxVar;
    mxVar(q2)=var(q2);
    if mod(i,10)==0
        imacs(mxCC);
        title(i);
        drawnow;
    end;
end;

%%
mxCC8=uint8(mxCC*100);
% save([basePath 'Info/' baseName 'ptmp.mat'],'mxCC8');
imwrite(rot90(mxCC8),[basePath 'Info/' baseName 'pcc.tif']);
return
%% Search for CC peaks

blankRadius=nt/2;
minAmp=1;
maxAmp=10;
amps=zeros(10,1);
mxCC2=mxCC;
nFound=0;
maxFound=50;
[amp ix iy]=max2di(mxCC2);
coords=[];
amps=zeros(2,1);
while amp>minAmp && nFound < maxFound
    % blank the cross-correlation
    if amp <= maxAmp && amp > minAmp  % valid peak
        nFound=nFound+1;
        coords(nFound,:)=[ix iy]; %-n/2-1;
        amps(nFound)=amp;
    end;
    mxCC2=mxCC2.*(1-fuzzymask(n,2,blankRadius,2,[ix iy]));
    imacs(mxCC2);
    title(nFound);
    drawnow;
    [amp ix iy]=max2di(mxCC2);
end;

[amps coords/2   ]


