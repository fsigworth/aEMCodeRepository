% GratingFresnel3.m
% Use the Fresnel propagator to model electron waves below the sample
backPropagationOn=1;
    trueLambda=.02;
    lambda=.25; % simulated wavelength
%     lambda=.02; % simulated wavelength

    nz=4096*3;
    dz=lambda/12; % number of A per z unit
%     total height will be equivalent of lambda/trueLambda * nz * dz
    effDz=dz*lambda/trueLambda;
    
    zBlock=20; % effective height of grating display, in A.
%     zBlock=.2;
    nx=512; % x dimension. Assumed to be even.
    d=5; % grating period, in A
    dx=d/12; % A per x unit
    halfCycles=3.75; % cycles per half-block
    halfCycles=5.75; % cycles per half-block
    nBlock=min(nx/2,halfCycles*d/dx); % half-width of block
%     nxx=nz;  % big square display
%     ctrx=floor((nxx+1)/2);
    
%     zBlock=200/dz*trueLambda/lambda; % no. angstroms thickness of specimen block

zFraction=.1; % grating distance from top.
if backPropagationOn
    zFraction=.3;
    zBlock=0;
end;
fpsi=zeros(nx,nz);
origZ=round(zFraction*nz);
origX=nx/2+1;
xiVals=(-nx/2:nx/2-1)';
ziVals=(nz-origZ:-1:-origZ+1)';

zCal=ziVals*effDz/10; % step in nm per z unit
xCal=xiVals*dx/10; % step in nm per x unit

amp=.4;
complexGrating=0;
weakPhase=0;
env=SquareWindow(nx,nx/16,1);
% env=Gaussian(nx,1,nx/4);
% env=ones(nx,1);
% envx=ones(nxx,1);

% Fill in unscattered wave
% zin=nx-origZ+1:nx;
% zVals=ziVals(nx-origZ+1:nx);
psi0=env*exp(1i*2*pi*dz*ziVals'/lambda); % undiffracted wave
% % psi0x=envx*exp(1i*2*pi*dz*ziVals'/lambda); % broad undiffracted wave


psi=psi0;


psi1=env;
block=(origX-nBlock:origX+nBlock)';
if complexGrating
    grating=exp(1i*2*pi*(block-origX)*dx/d);
else
    grating=cos(2*pi*(block-origX)*dx/d);
end;
fullGrating=zeros(nx,1);
fullGrating(block)=grating;

fullGrating2=zeros(nx,nx);
fullGrating2(block,block)=repmat(grating,1,2*nBlock+1);


if weakPhase
        psi1(block)=psi1(block).*(1+1i*amp*grating);
else
        psi1(block)=psi1(block).*exp(1i*amp*grating);
end;
% psi0=psi0.*Gaussian(nx,1,nx/4);
% psi0=psi0.*xWindow;
% psi(:,origZ)=1;
% psi(:,origZ)=1;
fpsi0=fftshift(fftn(ifftshift(psi1)));
fpsi(:,nz-origZ+1)=fpsi0;
psi(:,nz-origZ+1)=fftshift(ifftn(ifftshift(fpsi0)));

% 
%% propagate forward
for zi=1:nz-origZ % zi value
    zin=nz-origZ-zi+1; % z index in array, starts at nz-origZ, one below start.
    zStep=dz*zi;
%     zVal=zi*dz;
    H=exp(1i*2*pi*zStep/lambda)/(zStep).^0 ...
         .*exp(-1i*pi*lambda*zStep/(dx*nx)^2*(xiVals.^2));
%     
%       fpsi(:,zin)=(H.*fpsi(:,zin+1));
      fpsi(:,zin)=(H.*fpsi0);
    psi(:,zin)=ifftshift(ifftn(fftshift(fpsi(:,zin))));
end;

% Propagate backwards
if backPropagationOn
for zi=-1:-1:-origZ+1
    zin=nz-origZ-zi+1;
    zStep=dz*zi;
%     zVal=zi*dz;
    H=exp(1i*2*pi/lambda*zStep)/(zStep).^0 ...
        .*exp(-1i*pi*lambda*zStep/(dx*nx)^2*(xiVals.^2));
    
%       fpsi(:,zin)=(H.*fpsi(:,zin+1));
      fpsi(:,zin)=(H.*fpsi0);
    psi(:,zin)=ifftshift(ifftn(fftshift(fpsi(:,zin))));
end;    
end;
% 

%% Displays
% mysubplot(121);
% imacsx(fpsi,.3);
fs=18;
oneSlice=0;

figure(1);
set(gcf,'color',[.2 .2 .2]);
subplot(141);
% Mark the specimen position
amp=.7;
psid=psi;
nBlockZ=2*ceil(zBlock/(2*effDz));
zOrigi=nz-origZ;
blockZs=zOrigi+1:zOrigi+nBlockZ;
blockVals=repmat(fullGrating*.5+.5,1,nBlockZ);
psid(:,blockZs)=amp*blockVals+1i*inf;
% q1=psid(:,nz-origZ+1);
% psid(:,nz-origZ+1)=real(q1)/4+1i*imag(q1)/2;
pars=struct;
pars.x=xCal;
pars.x=xCal*10;
pars.y=[zCal(1) zCal(end)];
pars.sat=1;
pars.antiAliasing=1;
% imacx2(psid(:,zOrigi-1000:zOrigi+100),1,pars);
imacx2(psid,1,pars);
axis ij
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
ylabel('displacement z, nm');
xlabel('x, nm');
xlabel('x, angstroms');
title('\Psi(x,z)','color','w');

subplot(142);
psidNorm=psi./psi0;
psidNorm(:,zOrigi+1:zOrigi+nBlockZ)=amp*repmat(fullGrating*.5+.5,1,nBlockZ)+1i*inf;
pars.sat=1;
pars.scl=[];
pars.antiAliasing=0;
imacx2(psidNorm,.5,pars);
axis ij
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
ylabel('displacement z, nm');
xlabel('x, nm');
xlabel('x, angstroms');
title('\Psi/\Psi_0','color','w');

subplot(143);
disPsi=(psi-1*psi0)./psi0;
disPsi(:,blockZs)=.7*blockVals+1i*inf;
bkd=.5;
pars.antiAliasing=0;
pars.sat=1.5;
pars.scl=2;
if oneSlice
%     xctr=floor(nx/2+1);
%     xSlice=xctr+round(nBlock)+(-ceil(d/8):ceil(d/8))';
%     msk=ones(nx,nz,'single')*bkd;
%     msk(xSlice,:)=1;
    msk=1;
else
    msk=1;
end;
if complexGrating
    pars.scl=1.8;
end;
imacx2(disPsi.*msk,1,pars);
axis ij
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
ylabel('displacement z, nm');
xlabel('x, nm');
xlabel('x, angstroms');
if complexGrating
    title('\Psi_+','color','w');
else
    title('\Psi_+ + \Psi_-','color','w');
end;
pars.scl=[];
subplot(144);
psidr=psid;
q=isinf(imag(psidr));
psidr(q)=real(psidr(q))+.25;
% imags(pars.x/10,pars.y,abs(psidr))  % for nm X scale.
imags(pars.x,pars.y,abs(psidr))
axis ij
% imags(abs(psidr));
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
ylabel('displacement z, nm');
xlabel('x, nm');
xlabel('x, angstroms');
title('abs(\Psi)','color','w');
% figure(2);
% plot([real(q1) imag(q1) abs(q1)]);

return

%% Variable defocus movie

makeVideo=0;
vName='DefocusSim_wide_bprop_rev';

cFunction=2*sin(pi*zCal*10*trueLambda/d^2);
bkd=psidr*0+50;
xMag=xCal(end)*2;
nxd=min(nx,256);
% subplot(143);
% set(gca,'color',[.1 .1 .1]);
% hold on;
% plot(cFunction,zCal,'color',[1 1 1]);
% plot(zCal*0,zCal,'color',[0 0 1]);
% hold off;
% axis([-4 4 -inf inf]);

gMask=fuzzymask(nx,2,[1e4 nBlock*.8],[1e4 nBlock*.4]);
dGrating=Crop(fullGrating2.*gMask,nxd);
fGrating2=fftn(dGrating);
subplot(223);
    imaga(pars.x,pars.x,dGrating*100+128);
    axis equal off
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
   
if makeVideo
    v=VideoWriter([vName '.avi']);
            v.FrameRate=30;
    open(v);
    disp(['Making movie file: ' vName '.avi']);
end;


for def=90:-1:-90
% for def=0:1:250
% for def=250:-1:0
    c=ifftshift(CTF(nxd,dx,trueLambda,.001*def,0,0,0,0,0));
    filtImg=real(ifftn(fGrating2.*c));
    subplot(221);
    imaga(pars.x,pars.x,filtImg*100+128);
    axis equal off
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);

zval=zOrigi-def*10/effDz+zOrigi;

    
    subplot(144);
        imaga(pars.x,pars.y,abs(psidr));
        axis ij
    hold on;
plot(cFunction*xMag,zCal,'color',[.8 .8 .8],'linewidth',2);
plot(zCal*0,zCal,'color',[.6 .6 1]);
[val,ind]=min(abs(def-zCal));
plot(cFunction(ind)*xMag,def,'yo','markersize',20);
hold off;
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
title('sin(\piz\lambdas^2)','color','w');
    axis off;
    
subplot(143);
imags(pars.x,pars.y,abs(psidr))
axis ij
% imags(abs(psidr));
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
ylabel('displacement z, nm');
xlabel('x, nm');
xlabel('x, angstroms');
title('abs(\Psi)^2','color','w');
    
    hold on;
    plot([xCal(1) xCal(end)]*10,[def def],'y-','linewidth',2);
    hold off;
%     title(def,'color','w');
    drawnow;
    
        if makeVideo
    f=getframe(gcf);
    writeVideo(v,f);
    end;

    
end;

if makeVideo
    close(v);
end;




return






%% Try zoom-in of panel 1
subplot(141);
axis off;
mag=1;
while mag<20
    nz1=round(nz/mag);
    top=round((nz-origZ)+origZ/mag);
    bottom=top-nz1+1;
    ctr=nx/2+1;
    left=round(ctr-(nx/2/mag));
    right=round(ctr+((nx-2)/2/mag));
    psids=psid(left:right,bottom:top);
    imacx2(psids,1,[],1,3);
    axis off;
    drawnow;
    mag=mag*1.05;
end;
return

%% splice psi into psi0x
figure(2);
psi0x(xiVals+ctrx,:)=psi;
imacx2(psi0x);
