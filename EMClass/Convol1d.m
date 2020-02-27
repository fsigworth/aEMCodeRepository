% Convolution1D.m

fontSize=14;

bkd=[0 .14 .2]; % = 0 36 52 in Keynote.
set(groot,'defaultFigureColor',bkd);
set(groot,'defaultFigureInvertHardCopy','off');
set(groot,'defaultAxesColor',bkd);
set(groot,'defaultAxesXColor','w');
set(groot,'defaultAxesYColor','w');
set(groot,'defaultAxesFontSize',fontSize);
set(groot,'defaultLineLineWidth',2);
return







% IFTShah.m
% Show the FT of the shah function.

makeVideo=0;
% nameEnd='_square.jpg';
movieBaseName='Convol_1';

nt=200;
ntd=[25 50 100 200];  % initial nt display
nxd=[1 1 4];
frameRates=[4 10 15 30];
gw=  nt; % gaussian width
cd('/Users/fred/Documents/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.3bFT1D/figs/')
lw=2; % linewidth

% Make a dummy movie
% v=VideoWriter(movieBaseName);
%         open(v);

b1=60; % left/lower border
b2=30; % other border
bs=b1+b2;
nx=1920-bs;
nz=1080-bs;

figure(1);
clf;
h=gcf;
h.Color=bkd;

pos=h.Position;
h.Position=[pos(1:2) nx+bs nz+bs];


termColor=[.8 .8 .8];
tfs=18;  % text label fontsize
figure(1);
clf;
set(gcf,'color',bkd);
set(gcf,'InvertHardCopy','off');
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
set(gca,'color',bkd);
set(gca,'XColor','w');
set(gca,'YColor','w');
set(gca,'fontsize',18);
set(gca','YTickLabel',{});
ax2=gca;
subplot(212);
cla;
set(gca,'color',bkd);
set(gca,'XColor','w');
set(gca,'YColor','w');
set(gca,'fontsize',18);
ax1=gca;

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
            vName=[movieName num2str(itd)];
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
    hold on;
    plot(u,G,'w-','linewidth',lw);
    hold off;
    xlabel('u');
    ylabel('G(u)');
    
    axis([-uMax uMax 0 1]);
    
    
    subplot(212);
    cla;
    set(gca,'color',bkd);
    set(gca,'XColor','w');
    set(gca,'YColor','w');
    set(gca,'fontsize',18);
    ylabel('g(x)');
    ax1=gca;
    hold on;
    plot(x,y,'color',[1 1 1],'linewidth',lw);
%     plot(x,base,'--','color',[.3 .3 .3]);
    hold off;
    %         plot(x, y*0,'-','color',[.3 .3 .3],'linewidth',1.5);
    %         text(0.5,.7,num2str(val(i),3),'color','w','fontsize',18);
    %         for j=1:i
    %             text(x2s(1)+0.1,bases(1,j),num2str(2*fy(j+1),3),'verticalalignment','bottom',...
    %                 'color','w','fontsize',tfs);
    mx=max(20,max(y));
    axis([-inf inf -mx*.2 mx*1.2]);
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
print('ShahFTEnd.jpg','-djpeg');

