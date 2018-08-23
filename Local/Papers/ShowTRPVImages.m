% ShowTRPImages.m

cd('/Users/fred/EMWork/EMPIAR/10005/')
% cd('/Volumes/IndyT3/data/TRPV1')
% Load the first 1000 particles, sum of all frames
[mp,sp]=ReadMRC('particles/tv1_f01-30first5000.mrc',1,5000);
figure(1);
ImagicDisplay3(BinImage(mp,8));

%%
mc=ReadMRC('summed_frames/stack_0021_2x_SumCorr.mrc');
mc=RemoveOutliers(mc);
%
figure(2);
imags(GaussFilt(mc,.1));

%% Get CTF
Pa=struct;
Pa.defocus=1:.2:3;
Pa.deltadef=-.3:.1:.3;
Pa.theta=0:.1:pi/2;
Pa.lambda=EWavelength(300);
Pa.Cs=2;
Pa.B=10;
Pa.alpha=.1;
opts.blockSize=1024;
P=FitCTFk2(mc,Pa,sp.pixA,1,[.1 .3],opts);
% P=FitCTF(mc,Pa,sp.pixA,5,10);

%%


% Find particles in micrograph
mcIndex=22;
name=sprintf('stack_%04d_2x_SumCorr.mrc',mcIndex);

binFactor=4;  % downsampling to speed up search for particles

mc1=Crop(mc-mean(mc(:)),3840);
mcd=BinImage(mc1,binFactor);
fmcd=fftn(mcd);
np=50;
for i=1:np
    mp1=ifftshift(Crop(BinImage(mp(:,:,i),binFactor),3840/binFactor));
    cc=real(ifftn(fmcd.*conj(fftn(-mp1))));
    [val,ix(i),iy(i)]=max2d(cc);
end;
ix=ix*binFactor;
iy=iy*binFactor;
%%
imags(GaussFilt(mc1,.1))
hold on;
plot(ix,iy,'ys','markersize',32);
hold off;


%% Get map
[map,sm]=ReadMRC('EMD-5778/map/emd_5778.map');
map2=Downsample(map,128);

figure(3);
ShowSections(map2,[],45)


%% Make projections
betaStep=2;
gammaStep=2;
map2c=Crop(map2,80);
betas=0:betaStep:180;
nBetas=numel(betas);
gammas=0:gammaStep:90;
nGammas=numel(gammas);

allBetas=repmat(betas,nGammas,1);
allGammas=repmat(gammas',1,nBetas);
na=numel(allBetas);
angles=[zeros(na,1) allBetas(:) allGammas(:)];
projs=rsMakeTemplates(angles,map2c);
nProjs=size(angles,1);

%%

pprojs=Crop(projs,128,1);
refs=real(ifft2(fft2(pprojs).*repmat(abs(c),1,1,nProjs)));
refsv=reshape(refs,n^2,nProjs);
npm=44;
partsv=reshape(mf2(:,:,1:npm),n^2,npm);
ccs=refsv'*partsv;
ccs3=reshape(ccs,nGammas,nBetas,npm);

%%
figure(4);

for i=1:npm
ccn=ccs3(:,:,i);
ccnmx=max(ccn(:));
ccnm=ccn-ccnmx;
ccnm=max(ccnm,-100);
xcc=exp(ccnm);
imacs(ccnm)
% contourf(ccnm');
    title(i);
    pause
end;


%% Get particle metadata

load('particles/tv1_relion_data.mat');


%% load all the particles and downsample them
% mp0=ReadMRC('particles/tv1_f03-16.mrc');
% mp2=BinImage(mp0,2);
% mp2f=Downsample(mp0,128,1);
% WriteMRC(mp2f,s.pixA*2,'particles/tv1_f03-16ds2.mrc');
% clear mp0;

[mp2,s2]=ReadMRC('particles/tv1_f03-16ds2.mrc');

%% do phase flipping
mf2=mp2;
n=size(mp2);
nim=n(3);
n=n(1);
P.lambda=EWavelength(300);
mcIndex=0;
for i=1:nim
    ind=q.rlnMicrographName(i);
    if ind~=mcIndex
        mcIndex=ind;
        P.defocus=5e-5*(q.rlnDefocusU(ind)+q.rlnDefocusV(ind));
        P.deltaDef=5e-5*(q.rlnDefocusU(ind)-q.rlnDefocusV(ind));
        P.theta=q.rlnDefocusAngle(i)*pi/180;
        P.Cs=2;
        P.B=50;
        P.alpha=.01;
        c=-sign(ifftshift(CTF(n,s2.pixA,P)));
    end;

    mf2(:,:,i)=real(ifftn(fftn(mp2(:,:,i)).*c));
end;
%%  Sort particles by angles

tilt=10;
rot=45;

cls=abs(q.rlnAngleTilt-tilt)<6 & abs(q.rlnAngleRot-rot)<2;

nim=sum(cls)
xs=q.rlnOriginX(cls);
ys=q.rlnOriginY(cls);
clps=mf2(:,:,cls);
alphas=q.rlnAnglePsi(cls);
csh=clps;
for i=1:nim
    csh(:,:,i)=circshift(clps(:,:,i),round([xs(i) ys(i)]/2));
end;
csh=rsRotateImage(csh,-90-alphas);
figure(4);
subplot(221);
imags(rsMakeTemplates([180 tilt rot],map2));
subplot(223);
imags(sum(csh,3));


%% Make rotated, shifted, phase-flipped particle stack
nim=size(mf2,3);
ma2=mf2;
    xs=q.rlnOriginX;
ys=q.rlnOriginY;
alphas=q.rlnAnglePsi;
disp('shifting');
for i=1:nim
    ma2(:,:,i)=circshift(mf2(:,:,i),round([xs(i) ys(i)]/2));
end;
disp('rotating');
ma2=rsRotateImage(ma2,-90-alphas);  % in-plane aligned stack

 
    
    




%% Get corresponding projection for a particle

for ip=1:100
mp=mp2(:,:,ip);
xs=q.rlnOriginX(ip);
ys=q.rlnOriginY(ip);
mp=circshift(mp,round([xs ys]/2));
angs=[-90+q.rlnAnglePsi(ip) q.rlnAngleTilt(ip) q.rlnAngleRot(ip)];
proj=rsMakeTemplates(angs,map2);

        P.defocus=5e-5*(q.rlnDefocusU(ind)+q.rlnDefocusV(ind));
        P.deltaDef=5e-5*(q.rlnDefocusU(ind)-q.rlnDefocusV(ind));
        P.theta=q.rlnDefocusAngle(i)*pi/180;
        P.Cs=2;
        P.B=50;
        P.alpha=.1;
        c=-(ifftshift(CTF(n,s2.pixA,P)));
cproj=real(ifftn(fftn(proj).*c));

figure(4);
subplot(221);
imags(GaussFilt(mp,.2));
title(ip);
subplot(223);
imags(cproj);
drawnow;
pause
end;


%%  Sort particles by 3 angles

tilt=85;
rot=0;
psi=-90;

cls=abs(q.rlnAngleTilt-tilt)<6 & abs(q.rlnAngleRot-rot)<6 ...
        & abs(q.rlnAnglePsi-psi)<6;

nim=sum(cls)
xs=q.rlnOriginX(cls);
ys=q.rlnOriginY(cls);
clps=mf2(:,:,cls);
alphas=q.rlnAnglePsi(cls);
csh=clps;
for i=1:nim
    csh(:,:,i)=circshift(clps(:,:,i),round([xs(i) ys(i)]/2));
end;
% csh=rsRotateImage(csh,-90-alphas);
figure(4);
subplot(221);
imags(rsMakeTemplates([-90+psi tilt rot],map2));
subplot(223);
imags(sum(csh,3));

