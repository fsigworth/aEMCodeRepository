% GratingFresnelCxMovie.m
% Like the variable defocu movie in GratingFresnel3, but with complex CTF and image at
% diffraction plane

makeVideo=0;
vName='DefocusSim_complex_wide_bprop_rev';
epsi=0.1;

% Plot of the sin function of z
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

% Show the true object
gMask=fuzzymask(nx,2,[1e4 nBlock*.8],[1e4 nBlock*.4]);
dGrating=Crop(fullGrating2.*gMask,nxd);
fGrating2=fftn(ifftshift(dGrating));
% subplot(223);
    mysubplot(331);
    imaga(pars.x,pars.x,dGrating*100+128);
    axis equal off
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);
   
freqs=RadiusNorm(nxd)/dx;
CPars.lambda=trueLambda;
CPars.Cs=0;
CPars.B=0;
CPars.alpha=0


if makeVideo
    v=VideoWriter([vName '.avi']);
            v.FrameRate=30;
    open(v);
    disp(['Making movie file: ' vName '.avi']);
end;

% Defocus image
for def=90:-1:-90
% for def=0:1:250
% for def=250:-1:0
    c=ifftshift(CTF(nxd,dx,trueLambda,.001*def,0,0,0,0,0));
    filtImg=real(ifftn(fGrating2.*c));
    mysubplot(334);
    imaga(pars.x,pars.x,filtImg*100+128);
    axis equal off
set(gca,'fontsize',fs);
set(gca,'XColor',[1 1 1]);
set(gca,'YColor',[1 1 1]);

% Diffraction plane image
CPars.defocus=.001*def;
cx=ifftshift(ContrastTransfer(freqs,CPars,1));
diffImg=ifftshift(fGrating2).*(1+epsi*cx);
mysubplot(337);
imacx2(diffImg);
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


