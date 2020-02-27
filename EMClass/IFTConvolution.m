% IFTConvolution

figureScaling=.7;

ts=(-2:.01:2)';
us=ts*100;
a=10;
fs=18;

vlSetDarkGraphics(18);
co=get(groot,'DefaultAxesColorOrder');
co(1,:)=1;
set(groot,'DefaultAxesColorOrder',co);

vlSet1080Figure(1,figureScaling);

ys=ts*0;
ys(201)=100;
ys(51:150)=1;
ys(251:370)=randn(120,1);
ys=GaussFilt(ys,.15);

subplot(321);
plot(ts,ys);
axis([-inf inf -1.5 2]);
xlabel('x');
text(-1.5,1.6,'g(x)','fontsize',fs,'color','w');

subplot(322);
fy=fft(ifftshift(ys))/100;
plotcx(us,fftshift(fy));
xlabel('u');
text(-150,1.6,'G(u)','fontsize',fs,'color','w');

h=a*exp(-pi*(ts*a).^2);
subplot(323);
plot(ts,h);
xlabel('x');
ylabel('h');
text(-1.5,7,'h(x)','fontsize',fs,'color','w');

fh=fft(ifftshift(h))/100;
subplot(324)
plotcx(us,fftshift(fh));
xlabel('u');
text(-150,.75,'H(u)','fontsize',fs,'color','w');

fyFilt=fh.*fy;
yFilt=real(fftshift(ifft(fyFilt))*100);
subplot(325);
plot(ts,yFilt);
axis([-inf inf -1.5 2]);
xlabel('x');
text(-1.5,1.6,'g*h','fontsize',fs,'color','w');

subplot(326);
plotcx(us,fftshift(fyFilt));
xlabel('u');
text(-150,1.6,'G(u)H(u)','fontsize',fs,'color','w');

%% ------------------------ Deconvolution

vlSet1080Figure(2,figureScaling);

subplot(321);
plot(ts,[yFilt ys]);
cs=get(gca,'colororder');
plot(ts,ys,'color',cs(2,:));
hold on;
plot(ts,yFilt,'color','w');
hold off;
axis([-inf inf -1.5 2]);
xlabel('x');
text(-1.5,1.6,'Filtered g(x)','fontsize',fs,'color','w');

subplot(322);
plotcx(us,fftshift(fh.*fy));
xlabel('u');
text(-150,1.6,'Filtered G(u)','fontsize',fs,'color','w');

epsi=1e-7;

fhInv=1./(fh+epsi);
subplot(324);
semilogy(us,real(fftshift(fhInv)));
xlabel('u');
text(-150,.03/epsi,'H''(u)','fontsize',fs,'color','w');

hInv=fftshift(ifft(fhInv));
subplot(323);
plot(ts,hInv);
xlabel('x');
text(-1.5,3e6,'h''(x)','fontsize',fs,'color','w');

fyInv=fhInv.*fh.*fy;
subplot(326);
plotcx(us,fftshift(fyInv))
xlabel('u');
text(-150,1.6,'Deconvoluted G(u)','fontsize',fs,'color','w');


yFiltInv=real(fftshift(ifft(fyInv)));
subplot(325);
plot(ts,[yFiltInv*100 ys]);
xlabel('x');
text(-1.5,1.6,'Deconvoluted g(x)','fontsize',fs,'color','w');

% plot(ts,ys,'color',[1 .5 .5]);
% hold on;
% plot(ts,yFiltInv*100,'color','w');
% hold off;

axis([-inf inf -1.5 2]);

%% final illustration
nx=100;
ms=30;
fc=.1;
xs=(-1:1/nx:1);
nxc=round(1/(2*fc));

xss=1:nxc:2*nx+1;

vlSet1080Figure(3,figureScaling);

ys=SharpFilt(randn(2*nx+1,1),.1);
ys=(GaussFilt(ys,.1));
% yss=ys(1:1/(fc*2):2*nx+1);
% xss=xs(1:1/(fc*2):2*nx+1);
plot(xs,ys,'color',[.5 .5 1]);

hold on;
plot(xs(xss),ys(xss),'w.','markersize',ms);
hold off;

ax1=gca;
ax1.XTickLabels={};
ax1.YTickLabels={};


