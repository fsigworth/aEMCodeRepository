% MissingWedge.m

n=80;

pa='/Users/fred/EMWork/Simulations/GroEL/';
[map2,pixA]=ReadEMFile([pa 'GroelMap2A.mrc']);
n0=size(map2,1);
map4=Downsample(map2,n0/2);
map4c=Crop(map4,n);
rvol=map4c;

map2c=Crop(map2,2*n);
rvol=map2c;
n=2*n;


fvol=fftn(map2c);

% Create the missing wedge
[r,theta]=RadiusNorm(n);
thetaD=theta*180/pi;
sector=single((thetaD<-120)|(thetaD>-60 & thetaD<60) | (thetaD>120));
wedge=repmat(sector,1,1,n);
wedge=shiftdim(wedge,2);  %% make wedge X-Z
% ShowSections2(wedge);
sWedge=ifftshift(wedge);
wvol=real(ifftn(sWedge.*fvol));
ShowSections2(wvol);

halfPlane=ones(n,n,n);
halfPlane(:,:,1:n/2)=0;
outvol=wvol(:,n:-1:1,:).*mask3d.*halfPlane;
WriteMRC(outvol,pixA,[pa 'GroelMissingWedgeXZ.mrc']);





return;




% % OversampleTest
% Try OversamplingProject and OversamplingReconstruct functions

oversampleProjection=1;  % use the fast OversamplingProject function

ovs=4;  % Oversampling factor
nx=ovs*n;

sigmaN=2;   % noise added to each projection
useCTFs=0;  % give a random defocus to each projection


figure(1);
ShowSections2(rvol);
title('Original volume');
drawnow;

% angles=SphereAngles3(36,18,2);  % Angles in degrees
% angles=SphereAngles3(36,36,2);  % Angles in degrees
angles=(0:5:355)';
angles=(30:2:150)'
angles(:,2)=90;
angles(:,3)=0;

nim=size(angles,1);
disp(['making ' num2str(nim) ' projections']);

if oversampleProjection
    %     angs=[angles(:,3) angles(:,2) angles(:,1)];
    projs=OversamplingProject(rvol,angles);
    figure(6)
else  % Use the slower gridding function
    ks=3;
    vol=gridMakePaddedFT(rvol);
    comp=gridMakePreComp(n,ks);
    projs=(zeros(n,n,nim));
    angs=angles*pi/180;
    for i=1:nim
        p2=gridExtractPlane(vol,angs(i,:),ks);
        projs(:,:,i)=gridRecoverRealImage(p2,comp);
    end;
    figure(5);
end;

% Create the random ctfs.
if useCTFs
    ctfs=zeros([n n nim],'single');
    ctfz=zeros([n n nim],'single');
    defocus=1+rand(nim,1);
    % Simulate 4A pixel size, B=100
    for i=1:nim
        c=CTF(n,4,.025,defocus(i),2,100,.07);
        ctfz(:,:,i)=c;              % Zero frequency is in the center
        ctfs(:,:,i)=ifftshift(c);   % Zero frequency is at (1,1)
    end;
else
    ctfz=ones([n n nim],'single');
    ctfs=ctfz;
end;

projs=real(ifft2(fft2(projs).*ctfs))+sigmaN*randn(size(projs));

ImagicDisplay1(projs);
drawnow;

%%

% We do Fourier insertions into two volumes.  vol0 just has constant (or
% ctf^2) planes inserted, and is used for normalization.  vol1 is the actual
% reconstruction volume.  A Wiener filter then estimates vol1./vol0
sliceMask=fuzzymask(n,2,0.45*n,.05*n);  % Fourier mask for slices
slices0=ctfz.^2.*repmat(sliceMask,[1 1 nim]);
disp('First reconstruction');
% tic
dvol0=OversamplingReconstruct(slices0,angles);
% toc

figure(2);
ShowSections2(abs(dvol0).^.3);  % Show the Fourier volume
drawnow;

% Transforming projections
disp('Filtering the projections');
rmsk=fuzzymask(n,2,0.45*n,.05*n);  % real-space mask
slices=zeros([n n nim],'single');
for i=1:nim
    slices(:,:,i)=ctfz(:,:,i).*fftshift(fftn(ifftshift(rmsk.*projs(:,:,i))));
end;
disp('Second reconstruction');
tic
dvol1=OversamplingReconstruct(slices,angles);
toc
%%
% Wiener filter "normalization"
epsi=.0001;
fv0=fftn(dvol0);
dvol=fftn(dvol1).*fv0./(epsi+fv0.^2);

% Get the reconstructed volume.
disp('Final ifft');
revol=fftshift(real(ifftn(dvol)));
figure(3);

mask3d=fuzzymask(n,3,n*0.45,.1);
ShowSections2(revol.*mask3d);
title('Reconstructed volume');
drawnow;

halfPlane=ones(n,n,n);
halfPlane(n:-1:1,:,1:n/2)=0;
WriteMRC(revol.*mask3d.*halfPlane,pixA,[pa 'GroelMissingWedge.mrc']);

