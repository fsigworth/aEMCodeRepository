% TestFakeMicrograph
%
% Make a simulated micrograph with multiple vesicles.
%
deterministic=2;
useMask=0;

% Generic inverse filter params
f0=.02; f1=.0015; a1=.4; ex=4;

basePath='/Users/fred/matlabWork/Yunhui2/VesicleFinderData/120711/';
miName='Info/004_sq02_1_04mi.mat';
S0=.06;  % nominal scale
minS=.055;
maxS=.075;

load([basePath miName]);
imgName='Merged/004_sq02_1_04m.mrc';
ximg=ReadEMFile([basePath imgName]);

% Load the downsampled map and create a 2k x 2k image
load /Volumes/TetraData/Structures/AMPAR/3KG2map5.8AGoodScale.mat
nt=size(map,1);
map=map*5.8;  % V-A scaling
membraneOffset=-24/2;  % downsampled map by 2.

filtMap=map;
% filtMap=GaussHP(map,.05);  % 100 Å highpass

if deterministic==1
    mi.imageSize=[1024 1024];
    mi.vesicle.x=[256 512 768];  % zero-based values
    mi.vesicle.y=[512 512 512];
    mi.vesicle.r=[100 80 100];
    mi.vesicle.s=[.05 .05 .05];
    mi.vesicle.ok=[1 1 1];
    pAngles=[-90 90 0; 0 0 0; 90 90 90];
    isos=[0 0 0];
elseif deterministic==2
    baseName='Info/004_sq02_Fake';
    miName=[baseName 'mi.mat'];
    
    load([basePath miName]);
    mi.imageSize=[4096 4096];
    nx=16;
    ny=2;
%     r0=150/mi.pixA;
%     rstep=r0/2/nx;     % step from 150 to 300
    r0=64;
    rstep=0;     % step from 150 to 300
    xstep=mi.imageSize(1)/nx;
    ystep=mi.imageSize(2)/ny;
    xs=zeros(1,nx*ny);
    mi.vesicle.x=xs;
    mi.vesicle.y=xs;
    mi.vesicle.r=xs;
    mi.vesicle.s=xs;
    pAngles=zeros(nx*ny,3);
    isos=xs;
    alpha0=-5;
    alphaStep=20;
    betaStep=10;
    beta0=0;
    gammaStep=0;
    for i=1:ny
        for j=1:nx
            k=j+nx*(i-1);
           mi.vesicle.x(k)=(j-.5)*xstep;
           mi.vesicle.y(k)=(i-.5)*ystep;
           mi.vesicle.r(k)=r0+(nx-j)*rstep;
           mi.vesicle.s(k)=S0;
           mi.vesicle.ok(k)=1;
           beta=mod(beta0+(j-1)*betaStep,180);  % doesn't go above 180.
           pAngles(k,:)=[alpha0 beta (j-1)*gammaStep];
           isos(k)=~mod(j,2);
        end;
    end;    
end;

ds=2;
n=mi.imageSize(1)/ds;
nv=numel(mi.vesicle.x);
minR=32;
% minS=.047;
% maxS=.07;
nvTotal=numel(mi.vesicle.x);

mi.vesicle.ok=mi.vesicle.s>minS & mi.vesicle.s<maxS & mi.vesicle.r>20;
vesIndices=find(mi.vesicle.ok);
r0=200;
nRand=100;
mParts=zeros(n,n);
sumNp=0;
allLocs=[];
% vesIndices=45;
for i=vesIndices
    r=mi.vesicle.r(i);
    % emulate Poisson by taking a large number of trials with low
    % probability
    mu=r/r0;
    for iiso=0:1
        if deterministic
            iso=isos(i);
            np=iiso==iso;
        else
            np=sum(rand(nRand,1)<mu/nRand);  % poisson no. of particles
            iso=iiso;
        end;
        if np>0
            if deterministic
                alphas=pAngles(i,1);
                betas=pAngles(i,2);
                gammas=pAngles(i,3);
            else
                alphas=360*rand(np,1);
                betas=180/pi*acos(1-2*rand(np,1));
                gammas=180*rand(np,1);
            end;
                angles=[alphas betas gammas];
            if iso
%                 angles=[alphas+180 180-betas gammas];
                r=mi.vesicle.r(i)/ds+membraneOffset;
            else
                r=mi.vesicle.r(i)/ds-membraneOffset;
            end;
%             i
%             r
%             angles
            ctr=[mi.vesicle.x(i) mi.vesicle.y(i)]/ds+1;  % shift to 1-based coords
            [parts locs]=MakeRSCImage(n,filtMap,angles,r,ctr,iso);
            allLocs=[allLocs ; locs];
            mParts=mParts+parts*mi.vesicle.s(i);
        end;
        
        sumNp=sumNp+np;
    end; % for iso
end;
%%
sumNp
% allLocs
mi.vesicle.x=mi.vesicle.x(vesIndices);
mi.vesicle.y=mi.vesicle.y(vesIndices);
mi.vesicle.s=mi.vesicle.s(vesIndices);
mi.vesicle.r=mi.vesicle.r(vesIndices);
mi.vesicle.ok=mi.vesicle.s>0;
vesIndices=1:numel(vesIndices);
mi.particle.x=ds*(allLocs(:,1)-1);  % zero-based coordinates
mi.particle.y=ds*(allLocs(:,2)-1);


%%
% Get the CTF for filtering.
f=RadiusNorm(n)/(mi.pixA*ds);  % frequencies for evaluating the CTF
hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
hinv(f==0)=0;

hinv=1;
H=ifftshift(meGetEffectiveCTF(mi,n,ds).*hinv);
mParticles=-real(ifftn(fftn(mParts).*H));  % Filter with the ctf.
% mVesicles=meMakeModelVesicles(mi,n,vesIndices,0); % no ctf
% mVesicles=real(ifftn(fftn(mVesicles).*H));
mVesicles=meMakeModelVesicles(mi,n);
mVesicles=real(ifftn(fftn(mVesicles.*ifftshift(hinv))));
%%
nbin=1;
sigmaN=0;
if useMask
    msk=meGetMask(mi,n);
else
    msk=1;
end;
trueImage=msk.*GaussHP(Crop(ximg,n),.02);
simImage=msk.*(mVesicles+mParticles+sigmaN*randn(n,n));
figure(1); SetGrayscale;
imacs(BinImage(mParticles+sigmaN*randn(n,n),nbin));
figure(2); SetGrayscale;
imac(imscale(BinImage(simImage,nbin),256,1e-5));
% figure(3); SetGrayscale; imacs(BinImage(trueImage,nbin));