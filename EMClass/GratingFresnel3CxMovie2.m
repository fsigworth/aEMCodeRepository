% GratingCTF2D.m
% Use the complex CTF to model the image and diffraction plane
n=1024;
nz=128;
fpsi=zeros(n);
[xiVals, yiVals]=ndgrid((-n/2:n/2-1)');

dx=.5; % x or y per unit
lambda=.02;
d=10;
nblock=n/4;

weakPhase=0;


zCal=ziVals*dz*lambda/trueLambda/10; % physical units, in nm
xCal=xiVals(:,1)*dx/10;

% sTheta=lambda/d = .15. displacement of 250A in dz of 1600

backPropagationOn=0;


amp=.2;
% env=SquareWindow(nx,nx/16,1);
% env=Gaussian(nx,1,nx/4);
env=ones(n);
env=fuzzymask(n,2,0.4*n,.2*n);
% Fill in unscattered wave
% zin=nx-origZ+1:nx;
% zVals=ziVals(nx-origZ+1:nx);
psi0=reshape(env(:)*exp(1i*2*pi*dz*ziVals'/lambda),[n nz]); % undiffracted wave

psi=psi0;


psi1=env;

xs=(-n(1)/2:n(1)/2-1)';
% block=(origX-nblock:origX+nblock)';

fullGrating=repmat(cos(2*pi*xs/d),1,n(2));
% fullGrating=repmat(exp(1i*amp*cos(2*pi*xs/d)),1,n(2));
grating=fullGrating.*Crop(SquareWindow([2*nBlock+1 2*nBlock+1],ceil(nblock/2)),n);

if weakPhase
    psi1=env.*(1+1i*amp*grating);
else
    psi1=env.*exp(1i*amp*grating);
end;
% psi0=psi0.*Gaussian(nx,1,nx/4);
% psi0=psi0.*xWindow;
% psi(:,origZ)=1;
% psi(:,origZ)=1;
fpsi1=fftshift(fftn(ifftshift(psi1)));
fpsi(:,:,nz-origZ+1)=fpsi1;
psi(:,:,nz-origZ+1)=fftshift(ifftn(ifftshift(fpsi1)));
