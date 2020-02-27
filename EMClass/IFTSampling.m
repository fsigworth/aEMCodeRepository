% IFTSampling.m
% Demonstrate the sampling theorem
makeVideo=0;
% nameEnd='_square.jpg';
movieBaseName='ConvolFT';

nt=200;
ntd=[25 50 100 200];  % initial nt display
nxd=[1 1 4];
frameRates=[4 10 15 30];
gw=  nt; % gaussian width
cd('/Users/fred/Documents/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.3bFT1D/figs/')
% bkd=[0 .14 .2]; % = 0 36 52 in Keynote.
% lw=2; % linewidth
vlSetDarkGraphics(18);
vlSet1080Figure(1,.7);
co=get(groot,'DefaultAxesColorOrder');
co(1,:)=1;
set(groot,'DefaultAxesColorOrder',co);
fs=18;


xs=(-4:.01:4);
xShah=mod(xs,.25)==0;
sFreq=4;
% xShah=mod(xs,.67)==0;
% sFreq=1.5;
du=.01;
us=(-2:.01:2)';
usx=(-6:.01:6)';


g=exp(-pi*(xs-2).^2)+exp(-pi*(xs+2).^2);
G=2*exp(-pi*us.^2).*cos(4*pi*us);
% G=real(fftshift(fft(ifftshift(g))));
% G=Crop(G,numel(us));
subplot(321);
plot(xs,g);
xlabel('x');
ylabel('g(x)');

subplot(322);
co=get(gca,'colororder');
co(1,:)=1;
set(gca,'colororder',co);
% plot(G);
plot(usx,Crop(G,numel(usx)));
xlabel('u');
ylabel('G(u)');

subplot(323);
ax3=gca;
plot(xs,g.*xShah);
ax3.YTickLabel={};
xlabel('x');
txt=['g(x) III( ' num2str(sFreq) 'x )'];
text(0,.8,txt,'color','w','fontsize',fs,'horizontalalignment','center');

ng=numel(g)-1;
gsx=Crop(g.*xShah,ng*4+1);
G4=real(fftshift(fft(ifftshift(gsx))))/4;
G1=Crop(G4,numel(G));
subplot(324);
us1=us*3.14;
plot(us1,G1);
axis([-6 6 -2 2]);
xlabel('u');
txt=['G(u)*III( u/' num2str(sFreq) ' )'];
text(-2-(sFreq-4)*.5,1.5,txt,'color','w','fontsize',fs,'horizontalalignment','center');
G2=G1*(4/sFreq);
G2(abs(us1)>sFreq/2)=0;
subplot(326);
plot(us1,G2);
xlabel('u');
ylabel('G''(u)');
axis([-6 6 -2 2]);

g1=real(fftshift(ifft(ifftshift(G2))));
subplot(325);
xsc=Crop(xs,ng/8+1);
g1c=Crop(g1,ng/8+1);
g1cs=12.5*g1c;
plot(8*xsc,g1cs);
axis([-4 4 min(g1cs) 1]);
xlabel('x');
ylabel('g''(x)');

return


termColor=[.8 .8 .8];
tfs=18;  % text label fontsize
figure(1);
clf;
% set(gcf,'color',bkd);
% set(gcf,'InvertHardCopy','off');
np=8000;
np2=2*np+1;
nc=2.5;  % number of cycles
w=np/nc;

x=(-np:np)'/w;
uw=np/nt;
u=(-np:np)'/uw;
uCtr=np+1;
y=ones(np2,1);
G=zeros(np2,1);

G(uCtr)=1;
subplot(211);
cla;
ax1=gca;
% set(gca,'color',bkd);
ax1.YAxis.Color=bkd;
% set(gca,'YColor','w');
% set(gca,'fontsize',18);
ax1.YTickLabel={};
% ax2=gca;
subplot(212);
cla;
% set(gca,'color',bkd);
% set(gca,'XColor','w');
% set(gca,'YColor','w');
% set(gca,'fontsize',18);
ax2=gca;

base=zeros(np2,1);
uMax=0;  % start off showing only this many u points
itd=0;
% iVals=[1:ntd(1)];
% for j=2:numel(ntd)
%     iVals=[iVals ntd(j-1)+nxd(j):nxd(j):ntd(j)];
% end;

clear v;

for i=1:nt
    if i>uMax  % change movie segment
        itd=itd+1;
        uMax=ntd(itd);
        if makeVideo
            if exist('v','var')
                close(v);
                print(['ShahFT' num2str(itd) '.jpg'],'-djpeg');

            end;
            vName=[movieBaseName num2str(itd)];
            v=VideoWriter(vName);
            v.FrameRate=frameRates(itd);
            open(v);
            disp(['Making movie file: ' vName '.avi']);
        end;
    end;
    
    a=exp(-pi*(i/gw)^2);
    G(uCtr+i*uw)=a;
    G(uCtr-i*uw)=a;
    y=y+a*2*cos(2*pi*i*x);
    
%     if mod(i,nxd(itd))~=0
%         continue;
%     end;
    subplot(211);
%     hold on;
    plot(u,G,'w-','linewidth',lw);
if itd>3
    text(uMax*.5,.6,'G(u) = env(u) III(u)','fontsize',24,'color','w');
end;
%     hold off;
    xlabel('u','fontsize',24);
    ylabel('G(u)','fontsize',24);
%     ax1.YAxis.Visible='off';
    ax1.YAxis.TickLabels={};
    axis([-uMax uMax 0 1.1]);
    
    
    subplot(212);
    cla;
%     set(gca,'color',bkd);
%     set(gca,'XColor','w');
%     set(gca,'YColor','w');
%     set(gca,'fontsize',18);
%     ax1=gca;
%     hold on;
    plot(x,y,'color',[1 1 1],'linewidth',lw);
%     plot(x,base,'--','color',[.3 .3 .3]);
%     hold off;
    %         plot(x, y*0,'-','color',[.3 .3 .3],'linewidth',1.5);
    %         text(0.5,.7,num2str(val(i),3),'color','w','fontsize',18);
    %         for j=1:i
    %             text(x2s(1)+0.1,bases(1,j),num2str(2*fy(j+1),3),'verticalalignment','bottom',...
    %                 'color','w','fontsize',tfs);
    mx=max(20,max(y));
    axis([-inf inf -mx*.2 mx*1.2]);
    xlabel('x','fontsize',24);
    ylabel('g(x)','fontsize',24);
    drawnow
    
    if makeVideo
        f=getframe(gcf);
        writeVideo(v,f);
    end;
end;
if makeVideo
    close(v);
end;
title('');
print('ShahFTbEnd.jpg','-djpeg');

