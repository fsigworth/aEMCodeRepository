% TestFakeMicrograph2.m
% Make a simulated micrograph of AMPARs with multiple vesicles.
% Creates the images mParticles, mVesicles and mNoise.
%
deterministic=2;  % 2: dense array of particles
useMask=0;

% Generic inverse filter params
f0=.02; f1=.0015; a1=.4; ex=4;

basePath='/Volumes/TetraData/EMWork/Hideki/121210/Box_antagonist2_slot3/';
outPath='/Volumes/TetraData/EMWork/Hideki/121210/Simulation/';

miName='Info/sq02_10000mi.mat';
%  we will use this as the basis for our new mi file.

S0=.003;  % nominal scale for vesicles
minS=.002;
maxS=.005;

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


%  filtMap=GaussHP(map,.05);  % 100 Å highpass
 filtMap=GaussHP(map,.01*pixA);  % 100 Å highpass

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
%     baseName='Info/004_sq02_Fake';
%     miName=[baseName 'mi.mat'];
%     
%     load([basePath miName]);
%     mi.imageSize=[4096 4096];
    nx=16;
    ny=8;
%     r0=150/mi.pixA;
%     rstep=r0/2/nx;     % step from 150 to 300
    r0=88;
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
    betaStep=90/nx;
    beta0=0;
    gammaStep=180/ny;
    for i=1:ny
        for j=1:nx
            k=j+nx*(i-1);
           mi.vesicle.x(k)=(j-.5)*xstep;
           mi.vesicle.y(k)=(i-.5)*ystep;
           mi.vesicle.r(k)=r0+(nx-j)*rstep;
           mi.vesicle.s(k)=S0;
           mi.vesicle.ok(k,:)=true(1,4);
           beta=mod(beta0+(j-1)*betaStep,180);  % doesn't go above 180.
           pAngles(k,:)=[alpha0 beta floor((i-1)/2)*gammaStep];
           isos(k)=~mod(i,2);
        end;
    end;    
end;

ds=2;
% n=mi.imageSize(1)/ds;
nv=numel(mi.vesicle.x);
minR=32;
% minS=.047;
% maxS=.07;
nvTotal=numel(mi.vesicle.x);

% mi.vesicle.ok=mi.vesicle.s>minS & mi.vesicle.s<maxS & mi.vesicle.r>20;
vesIndices=find(mi.vesicle.ok(:,2));
r0=200;
nRand=100;
mParts=zeros(n);
sumNp=0;
allLocs=[];
locs=single(zeros(1,8));
% vesIndices=45;
for i=vesIndices(:)'
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
            [parts, locs(1,1:2)]=MakeRSCImage(n,filtMap,angles,r,ctr,iso);
            locs(:,4)=i;  % mark the vesicle number
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
mi.vesicle.ok=true(numel(mi.vesicle.s),4);
mi.particle.picks=single(zeros(size(allLocs,1),8));
mi.particle.picks=ds*(allLocs-1);  % zero-based coordinates
mi.particle.picks(:,3)=16;  % set the manual pick flag

%%
% Get the CTF for filtering.
f=RadiusNorm(n)/(mi.pixA*ds);  % frequencies for evaluating the CTF
hinv=a1./(1+(f1./f).^2)+(1-a1)*1./(1+(f0./f).^ex);
hinv(f==0)=0;

% hinv=1;
H=ifftshift(meGetEffectiveCTF(mi,n,ds).*hinv);
mParticles=-real(ifftn(fftn(mParts).*H));  % Filter with the ctf.

mVesicles0=meMakeModelVesicles(mi,n);
% mVesicles0=meMakeModelVesicles(mi,n,0,0);  % no CTF version
mVesicles=real(ifftn(fftn(mVesicles0.*ifftshift(hinv))));
mi0=mi;
%%  Final images and display
mi=mi0;
nbin=2;
sigmaN=.3;  % noisy image
% sigmaN=.1;
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
figure(2); SetGrayscale;
imac(imscale(BinImage(simImage,nbin),256,1e-5));
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
