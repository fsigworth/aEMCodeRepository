function [imgs,si]=reMakeSimStackSub(sm)
% The simulation parameters are
% sm.mDefoci   Defocus for each micrograph, nMics x 1.
% sm.miIndex   Micrograph number for each image
% sm.angles    nImgs x 3
% sm.rocks     nImgs x 2, alpha,beta rocking
% sm.bobs      nImgs x 1
% sm.clicks    nImgs x 2
% sm.isos      nImgs x 2, boolean
% sm.vesR      nImgs x 1
% sm.pixA
% sm.mbnOffset 1x1, in pixels; zero to use map default.
% sm.sigmaN
% sm.imgAmp
% sm.n
% sm.mapName

n=sm.n;
pixA=sm.pixA;

% Working image sizes
nc=NextNiceNumber(1.2*n);  % padded image for CTF etc.
nMicrographs=numel(sm.mDefoci);

B0=50;
B1=50;
lambda=.025;

% Angle and shift search

nVols=1;  % only implemented 1 volume at present

[vols, mbnOffsetA]=arGetRefVolumes(pixA,n,sm.mapName,nVols);
if sm.mbnOffset==0   
    sm.mbnOffset=mbnOffsetA/pixA;
end;

nImgs=size(sm.angles,1);
si=struct;
si.miIndex=sm.miIndex;
si.miParticle=zeros(nImgs,1,'uint16');
si.alpha0=zeros(nImgs,1);
si.rVesicle=sm.vesR;
si.pixA=pixA;
si.mbnOffset=sm.mbnOffset;
si.weights=[1 1];
si.ctfs=zeros(n,n,nMicrographs,'single');

ctfsc=zeros(nc,nc,nMicrographs,'single');
for j=1:nMicrographs
    d=sm.mDefoci(j);
    si.ctfs(:,:,j)=abs(CTF(n,pixA,lambda,d,2,B0+B1*d,0));
    ctfsc(:,:,j)=abs(CTF(nc,pixA,lambda,d,2,B0+B1*d,0));
end;

% Make the images -------------
[imgs0,si.yClick,sm.partAngles]=rsMakeFakeImages2(vols,nc,si.rVesicle,sm);

msk=fuzzymask(nc,2,0.45*n,0.1*n);
imgs1=zeros(nc,nc,nImgs,'single');
for i=1:nImgs
    c=ctfsc(:,:,si.miIndex(i));
    imgs1(:,:,i)=msk.*real(ifftn(fftn(imgs0(:,:,i)).*ifftshift(c)));
end;

imgs=sm.imgAmp*Crop(imgs1,n,1)+sm.sigmaN*randn(n,n,nImgs);
