% TestFakeMicrographKv.m
% Make a simulated micrograph of Kv channels with multiple vesicles.
% Creates the images mParticles, mVesicles and mNoise.
% Based on TestFakeMicrograph5 is simplified, assumes the "deterministic=3" option of the previous version.
%
useMask=0;
freezeNoise=0;  % force re-use of the noise field

% % Generic inverse filter params
% f0=.02; f1=.0015; a1=.4; ex=4;

% basePath='/Volumes/TetraData/EMWork/Hideki/121210/Simulation/SimDefocusSeries/';
% outPath=basePath;
ds=4;  % downsample by this factor

defoci=[2 10];
doses=[15 20];
nZeros=1;

%  we will use this as the basis for our new mi file.

S0=.004;  % nominal scale for vesicles; conversion V� to contrast

    [nm, pathName]=uigetfile({'*mi.mat' '*mi.txt'},'Select mi file');
    if isnumeric(pathName) % File selection was cancelled
        return
    end;
    [rootPath, infoPath]=ParsePath(pathName);
    cd(rootPath);

    
mi=ReadMiFile([infoPath nm]);

n=mi.imageSize/ds;
ds=mi.imageSize(1)/n(1);
pixA=mi.pixA*ds;

if ~(freezeNoise && exist('mNoise','var'))  % Use a consistent noise field
    baseNoise=randn(n);
end;


mi.doses=doses;
nd=numel(defoci);
for i=1:nd
    mi.ctf(i)=mi.ctf(1);
    mi.ctf(i).defocus=defoci(i);
    mi.ctf(i).B=150+100*defoci(i);
end;
mi.ctf(nd+1:end)=[];
% mi.doses=10*ones(1,nd);
mi.weights=ones(1,nd);


% Get the AMPAR map
matName='~/aEMCodeRepository/Data/3KG2mapsub58.mat';

s=load(matName);
s.pixA=5.8;
origMap=s.sim;
membraneOffset=s.mbnOffsetA/pixA;
% [origMap, s]=ReadMRC(mapName);

% Fix the voxel size
mpixA=s.pixA;  % original voxel isze
nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
nt=ceil(nt1/8)*8;
[map, finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
map=map*pixA/mpixA;  % approx amplitude correction (V-A scaling)
membraneOffset=-70/pixA;  % downsampled map by 2.

% filtMap=GaussHP(map,.01*pixA);  % 100 � highpass
filtMap=map;

nx=10;
ny=8;
% nx=8;
% ny=6;
vesRadius=150;
%     r0=150/mi.pixA;
%     rstep=r0/2/nx;     % step from 150 to 300
r0=vesRadius/mi.pixA;  % radius in pixels
rstep=5;     % step from 150 to 300
xstep=mi.imageSize(1)/nx;
ystep=mi.imageSize(2)/ny;
xs=zeros(nx*ny,1);
mi.vesicle.x=xs;
mi.vesicle.y=xs;
mi.vesicle.r=xs;
mi.vesicle.s=xs;
mi.vesicle.ok=false(nx*ny,4);
pAngles=zeros(nx*ny,3);
isos=xs;
alpha0=0;
alphaStep=0;
beta0=90;
% betaStep=180/(ny+1);
betaStep=0;
betas=beta0+betaStep*((1:ny)-.5);
% gammaStep=180/(nx+1);
gamma0=0;
gammaStep=0;
gammas=gamma0+gammaStep*((1:nx)-.5);
for i=1:ny
    for j=1:nx
        k=j+nx*(i-1);
        mi.vesicle.x(k)=(j-.5)*xstep;
        mi.vesicle.y(k)=(i-.5)*ystep;
        mi.vesicle.r(k)=r0+(nx-j)*rstep;
        mi.vesicle.s(k)=S0;
        mi.vesicle.ok(k,:)=true(1,4);
        pAngles(k,:)=[alpha0 betas(i) gammas(j)];
        isos(k)=0;
    end;
end;

% ds=2;
% n=mi.imageSize(1)/ds;
nv=numel(mi.vesicle.x);
nvTotal=numel(mi.vesicle.x);

vesIndices=find(mi.vesicle.ok(:,2));
mParts=zeros(n);
sumNp=0;
allLocs=[];
locs=single(zeros(1,8));
% vesIndices=45;
for i=vesIndices(:)'
    r=mi.vesicle.r(i);
    iso=isos(i);
    alphas=pAngles(i,1);
    betas=pAngles(i,2);
    gammas=pAngles(i,3);
    angles=[alphas betas gammas];
    if iso
        r=mi.vesicle.r(i)/ds+membraneOffset;
    else
        r=mi.vesicle.r(i)/ds-membraneOffset;
    end;
    ctr=[mi.vesicle.x(i) mi.vesicle.y(i)]/ds+1;  % shift to 1-based coords
    [parts, locs(1,1:2)]=MakeRSCImage(n,filtMap,angles,r,ctr,iso);
    locs(:,4)=i;  % mark the vesicle number
    allLocs=[allLocs ; locs];
    mParts=mParts+parts*mi.vesicle.s(i);  % fractional contrast
end;

%%
% Mysterious scale factors to match merged image amplitudes.
particleScale=.25;
noiseScale=.15;

mi.particle.picks=single(zeros(size(allLocs,1),8));
mi.particle.picks=ds*(allLocs-1);  % zero-based coordinates
mi.particle.picks(:,3)=16;  % set the manual pick flag

mi.particle.true=mi.particle.picks;
mi.trueVesicle=mi.vesicle;
%
% Get the CTF for filtering.
% f=RadiusNorm(n)/(mi.pixA*ds);  % frequencies for evaluating the CTF
% hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
% hinv(f==0)=0;
% 
hinv=meGetNoiseWhiteningFilter(mi,n,ds);
% hinv=1; %------------ No noise whitening ***********-------------

H=ifftshift(meGetEffectiveCTF(mi,n,ds,0,nZeros).*hinv);
figure(3);
plot(H(1:n(1)/2,1));
mParticles=-real(ifftn(fftn(mParts).*H))*particleScale*doses(1);  % Filter with the ctf.

% Make the filtered vesicle image, scaled as fractional contrast
mVesicles0=meMakeModelVesicles(mi,n,0,0);  % no CTF version, fractional contrast
mVesicles=real(ifftn(fftn(mVesicles0).*H)); % filter with CTF and prewhitening
mi0=mi;
%  Final images and display
mi=mi0;
nbin=1;
fc=.25;

sigmaN=sqrt(doses(1))/pixA * noiseScale; % effective 'shot' noise, in fractional contrast

if useMask
    msk=meGetMask(mi,n);
else
    msk=1;
end;
% trueImage=msk.*GaussHP(Crop(ximg,n),.02);
    mNoise=sigmaN*baseNoise;
%     mNoise(1:n(1)/4,:)=0;  % erase noise from half the image
simImage=msk.*(mVesicles+mParticles+mNoise);
partImage=msk.*(mParticles+mNoise);

partImageFilt=GaussFilt(mParticles+mNoise,fc);
xBetas=(0:n(1)-1)*(180+betaStep)/n(1)-betaStep/2;
yGammas=(0:n(2)-1)*(180+gammaStep)/n(2)-gammaStep/2;

figure(1); SetGrayscale;
imac(xBetas,yGammas,imscale(BinImage(partImageFilt,nbin),256,1e-5));

figure(2); SetGrayscale;
% simImageFilt=GaussFilt(simImage,fc);
% % imac(xBetas,yGammas,imscale(BinImage(simImage,nbin),256,1e-5));
% imac(xBetas,yGammas,imscale(BinImage(simImageFilt,nbin),256,1e-5));
% xlabel('Beta');
% ylabel('Gamma');
imags(simImage);

figure(3);
imags(partImage);

% mi.basePath=outPath;
% mi.imagePath='Micrograph/';
% mi.procPath='Merged/';
% mi.infoPath='Info/';
% mi.baseFilename=[mi0.baseFilename 'Sim_' num2str(100*sigmaN)];
% outMiName=[outPath mi.infoPath mi.baseFilename 'mi.mat']
% save(outMiName,'mi');
% outMergeName=[outPath mi.procPath mi.baseFilename 'm.mrc']
% WriteMRC(simImage,pixA,outMergeName);
