% TestFakeMicrograph3.m
% Make a simulated micrograph of AMPARs with multiple vesicles.
% Creates the images mParticles, mVesicles and mNoise.
% version 3 is simplified, assumes the "deterministic=3" option of the previous version.
%
useMask=0;

% Generic inverse filter params
f0=.02; f1=.0015; a1=.4; ex=4;

basePath='/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/';
outPath='/Volumes/TetraData/EMWork/Hideki/121210/Simulation/';

miName='Info/sq02_10000mi.mat';
%  we will use this as the basis for our new mi file.

S0=.003;  % nominal scale for vesicles

load([basePath miName]);
ximg=meReadMergedImage(mi);
n=size(ximg);
ds=mi.imageSize(1)/n(1);
pixA=mi.pixA*ds;

mapName='/Volumes/TetraData/Structures/AMPAR/3KG2map58.mrc';

[origMap s]=ReadMRC(mapName);
mpixA=s.pixA;
nt1=size(origMap,1)*mpixA/pixA;  % final effective map size
nt=ceil(nt1/8)*8;
[map, finalmag]=DownsampleGeneral(origMap,nt,mpixA/pixA);
map=map*pixA;  % approx amplitude correction (V-A scaling)
membraneOffset=-70/pixA;  % downsampled map by 2.


filtMap=GaussHP(map,.01*pixA);  % 100 Å highpass
filtMap=map;

nx=15;
ny=12;
%     r0=150/mi.pixA;
%     rstep=r0/2/nx;     % step from 150 to 300
r0=88;  % radius in pixels
rstep=0;     % step from 150 to 300
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
betaStep=180/(ny+1);
betas=betaStep*((1:ny)-.5);
gammaStep=180/(nx+1);
gammas=gammaStep*((1:nx)-.5);
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

ds=2;
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
    mParts=mParts+parts*mi.vesicle.s(i);
end;

%%
mi.particle.picks=single(zeros(size(allLocs,1),8));
mi.particle.picks=ds*(allLocs-1);  % zero-based coordinates
mi.particle.picks(:,3)=16;  % set the manual pick flag

mi.particle.true=mi.particle.picks;
mi.trueVesicle=mi.vesicle;
%%
% Get the CTF for filtering.
f=RadiusNorm(n)/(mi.pixA*ds);  % frequencies for evaluating the CTF
hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
hinv(f==0)=0;

%         hinv=meGetNoiseWhiteningFilter(mi,n,ds);

H=ifftshift(meGetEffectiveCTF(mi,n,ds).*hinv);
mParticles=-real(ifftn(fftn(mParts).*H));  % Filter with the ctf.

mVesicles0=meMakeModelVesicles(mi,n);
% mVesicles0=meMakeModelVesicles(mi,n,0,0);  % no CTF version
mVesicles=real(ifftn(fftn(mVesicles0.*ifftshift(hinv))));
mi0=mi;
%%  Final images and display
mi=mi0;
nbin=2;
sigmaN=0;  % noisy image is 0.4

if useMask
    msk=meGetMask(mi,n);
else
    msk=1;
end;
trueImage=msk.*GaussHP(Crop(ximg,n),.02);
mNoise=sigmaN*randn(n);
simImage=msk.*(mVesicles+mParticles+mNoise);
figure(1); SetGrayscale;
imacs(BinImage(mParticles+mNoise,nbin));
xBetas=(0:n(1)-1)*(180+betaStep)/n(1)-betaStep/2;
yGammas=(0:n(2)-1)*(180+gammaStep)/n(2)-gammaStep/2;
figure(2); SetGrayscale;
imac(xBetas,yGammas,imscale(BinImage(simImage,nbin),256,1e-5));
xlabel('Beta');
ylabel('Gamma');
% figure(3); SetGrayscale; imacs(BinImage(trueImage,nbin));

mi.basePath=outPath;
mi.imagePath='Micrograph/';
mi.procPath='Merged/';
mi.infoPath='Info/';
mi.baseFilename=[mi0.baseFilename 'Sim_' num2str(100*sigmaN)];
outMiName=[outPath mi.infoPath mi.baseFilename 'mi.mat']
save(outMiName,'mi');
outMergeName=[outPath mi.procPath mi.baseFilename 'm.mrc']
WriteMRC(simImage,pixA,outMergeName);
