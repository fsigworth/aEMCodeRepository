% GratingFresnel2D.m
% Use the Fresnel propagator to model electron waves below the sample
n=256*[1 1]; % dimensions of outputs. Assumed to be even.
nz=1024;
zFraction=.2;
fpsi=zeros([n nz]);
origZ=round(zFraction*nz);
origX=n(1)/2+1;
[xiVals yiVals]=ndgrid((-n/2:n/2-1)');
ziVals=(nz-origZ:-1:-origZ+1)';

dx=.5; % x or y per unit
dz=1/12; % z per unit
lambda=1; % 10x actual
trueLambda=.02;
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

% 
%% propagate forward
for zi=1:nz-origZ % zi value
    zin=nz-origZ-zi+1; % z index in array, starts at nz-origZ, one below start.
    zStep=dz*zi;
%     zVal=zi*dz;
    H=exp(1i*2*pi*zStep/lambda) ...
         .*exp(-1i*pi*lambda*zStep/(dx*n(1))^2*(xiVals.^2+yiVals.^2));
%     
%       fpsi(:,zin)=(H.*fpsi(:,zin+1));
      fpsi(:,:,zin)=(H.*fpsi1);
    psi(:,:,zin)=ifftshift(ifftn(fftshift(fpsi(:,:,zin))));
end;

% Propagate backwards
if backPropagationOn
for zi=-1:-1:-origZ+1
    zin=nz-origZ-zi+1;
    zStep=dz*zi;
%     zVal=zi*dz;
    H=exp(1i*2*pi/lambda*zStep)/(zStep).^0 ...
        .*exp(-1i*pi*lambda*zStep/(dx*n)^2*(xiVals.^2+yiVals.^2));
    
%       fpsi(:,zin)=(H.*fpsi(:,zin+1));
      fpsi(:,:,zin)=(H.*fpsi1);
    psi(:,:,zin)=ifftshift(ifftn(fftshift(fpsi(:,:,zin))));
end;    
end;
% 

%%
% mysubplot(121);
% imacsx(fpsi,.3);

figure(1);

% Mark the specimen position
psid=squeeze(psi(:,origX,:));
psis=squeeze(psi0(:,origX,:));
q1=psid(:,nz-origZ+1);
psid(:,nz-origZ+1)=real(q1)/2+1i*imag(q1)/2;

subplot(141);
imacx2(psid,.5);

subplot(142);
imacx2((psid-psis)./psis,.5);

subplot(143);
imags(xCal,-zCal,abs(psid))

%%
figure(2);
close(2);
scale=2/3;
vlSetDarkGraphics(round(18*scale));
vlSet1080Figure(2,2/3);
% figure(2);
% clf;

msk=fuzzymask(256,2,90,40);
izStep=16;
mulr=64/amp;
addr=64;
offs=.9;
lz=140;
aMax=2800;
for iz=nz:-izStep:lz
subplot(121)
    imaga(mulr*(abs(psi(:,:,iz)./psi0(:,:,iz))-1)+addr);
subplot(122);
    dpsi=psi(:,:,iz)./psi0(:,:,iz)-offs;
    fpsi=fftshift(fftn(ifftshift(msk.*dpsi)));
% % %     attempt at a raised-cosine 
% %     a=abs(fpsi);
% %     s1=.1;
% %     mul=.5-.5*cos(s1*a);
% %     mul(a*s1>pi)=1;
%     aMax=max(a(:));
    mul=1;
    fpsi=fpsi.*mul;
    pars.scl=1.5/aMax;
     imacx2(Crop(fpsi,3*n/8),1,pars);
%     imacx2(dpsi);
    title(nz-iz);
    drawnow;
end;

% figure(2);
% plot([real(q1) imag(q1) abs(q1)]);
% 

