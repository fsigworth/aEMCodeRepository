% Orientability.m


saveProjections=0;

ICDOnly=1;  % align only the intracellular domain of Kv1.2
n=96;

load('/Users/fred/aEMCodeRepository/AMPAR/KvMap.mat');
pixA0=pixA;  % original map voxel, =2

pixA=2.494;  % size we will use
mapx=Crop(map,120);

if ICDOnly
    mapx(:,:,71:end)=0;
end;

offsets=CenterOfMass(mapx);
offsets

map1=DownsampleGeneral(circshift(mapx,[0 0 -round(offsets(3))]),n,pixA0/pixA);
ShowSections(map1,[],45);

%%
angles=SphereAngles3(9,18,1);  % Number of steps in psi, theta: 10 degrees, one hemisphere
angles(:,1)=angles(:,1)/4;  % psi reduced to 90 degree range.
% No in-plane rotations needed for the analysis.
% alpha=0:10:170;
% nAlpha=18;
dA=1;  % angle perturbation

n=size(map1,1);
nim=size(angles,1);

pAngles=zeros(nim,3,4);
pAngles(:,:,1)=angles;
pAngles(:,:,2)=angles+[dA 0 0];
pAngles(:,:,3)=angles+[0 dA 0];
pAngles(:,:,4)=angles+[0 0 dA];
disp(['making ' num2str(nim) ' projection triplets']);

ks=5;
vol=gridMakePaddedFT(map1);
comp=gridMakePreComp(n,ks);
projs=(zeros(n,n,nim,3,'single'));
pAngs=pAngles*pi/180;
for i=1:nim
    for j=1:4
        p2=gridExtractPlane(vol,pAngs(i,:,j),ks);
        projs(:,:,i,j)=gridRecoverRealImage(p2,comp);
    end;
end;
%
fps=fft2(projs);
diffs=zeros(n,n,nim,3);
sOrient=zeros(n/2,3);
sDiffs=zeros(n,n,3);

fs=(0:n/2-1)/(n*pixA);

%%
for j=1:3
    diffs(:,:,:,j)=projs(:,:,:,j+1)-projs(:,:,:,1);
    fDiffs=abs(fft2(diffs(:,:,:,j))).^2/n^2;
    sDiffs(:,:,j)=fftshift(mean(fDiffs,3));
    sOrient(:,j)=Radial(sDiffs(:,:,j));
end;

if saveProjections
    save Orient.mat fs sDiffs sOrient map1 pixA dA
    return
end;
%%

% load Orient.mat

% Load a stack, to get ctf parameters from one of the mi structures.
n=size(map1,1);
disp('Getting a stack file...');
[stackName,stackPath]=uigetfile('*si.mat');
disp(stackName);
[rootPath,stackDir]=ParsePath(stackPath);
cd(rootPath);
load([stackDir stackName]);

%%
ind=1;  % mi index to use
mi=si.mi{ind};
mi.weights(2)=1;
% mi.ctf(1).defocus=defoci(idef);
mi.ctf(1).B=40;
fs=(0:n/2-1)'/(pixA*n);
c2=sectr(meGetEffectiveCTF(mi,n,2)).^2;  % downsampling by 2
a3=8;
a2=1;
pars=[1 .31 6];
bkdSpectrum=1+a3*c2.*(CarbonSpectrum(fs,pars)-1)+a2*c2;
c2n=c2./bkdSpectrum;
%
%         c2=ContrastTransfer(fs,.025,defocus,2,B,.02).^2;
figure(3);
subplot(221)
csO=sOrient.*repmat(c2n,1,3);
plot(fs,csO);
css=sum(csO);
% csd(:,idef)=css(:);
legend(['  phi ' num2str(css(1),3)],...
    ['theta ' num2str(css(2),3)],...
    ['  psi ' num2str(css(3),3)]);
ylabel('Orientation density');
xlabel('Frequency A^{-1}');
[root,main]=ParsePath(mi.basePath);
title([main '  ' mi.baseFilename],'interpreter','none');
%     legend(num2str(sum(csO)',3));
subplot(223);
plot(fs,[c2 bkdSpectrum/10]);
legend('CTF^2','noise Spectrum');
xlabel('Frequency A^{-1}');

subplot(224);
plot(fs,c2n);
ylabel('CTF^2/noise');
xlabel('Frequency A^{-1}');

subplot(222);
plot(fs,cumsum(csO(:,2)));
ylabel('Theta orientation power');
xlabel('Frequency A^{-1}');
title(['weights [' num2str(mi.weights) ']  doses ' num2str(mi.doses,3)]);
% defs=0;
% for i=1:numel(mi.ctf)
%     defs(i)=mi.ctf(i).defocus;
% end;
% title(['Defocus   ' num2str(defs,3)]);
% subplot(2,2,2);
% imags(sDiffs(:,:,2));

drawnow;
