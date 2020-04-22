% IFTGaussianToDelta.m
% Make a movie of limiting Gaussians
% The movie is written out as a series of jpeg images.
% fs Nov 2018
writeImgs=0;
% nameEnd='_square.jpg';
figNameExt='_gauss.jpg';

vlSetDarkGraphics;
% bkd=[0 .14 .2]; % = 0 36 52 in Keynote.
% termColor=[.8 .8 .8];
brighter=0.3;
lw=2;
lw2=4;
tfs=18;  % text label fontsize
% fs=18;
figure(1);
clf;
set(gcf,'color',bkd);
% set(gcf,'InvertHardCopy','off');
np=200;
w=np/4;

xSpan=3;
xShift=1.5;
bShift=2;
fSpan=3;
colorOrderIndex=1;

xs=xSpan*(-2*w:2*w-1)'/w+xShift;
xs2=xs/2;
fs=fSpan*(-2*w:2*w-1)'/w;
fs2=fs/2;

ax1=subplot(321); % g(x)
% ax1.ColorOrderIndex=colorOrderIndex;
origColorOrder=ax1.ColorOrder;

% ax2=subplot(322); % for the FTs

bShift=2;

y=exp(-pi*(xs2-bShift).^2);

fAng=exp(-1i*2*pi*bShift*fs2);
fMag=exp(-pi*(fs2).^2);
fy=fMag.*fAng;


% aVals=[1 2 4 8 16];
% nt=numel(aVals);
% terms=zeros(numel(x2s),nt);
% iTerms=zeros(numel(x2s),nt);
yLims1=[0 1.3];
yLims2=[-1.3 1.3];

% Plot the original function
ax1=subplot(321);
plot(xs2, y,'linewidth',lw);
axis([xs2(1) xs2(end) yLims1]);
ax1.ColorOrder=min(origColorOrder+brighter,1);
ylabel('g(x)');
xlabel('x');



ax2=subplot(322); % plot the FTs

plot3(fs2,0*real(fy),0*imag(fy),'-','color',[.7 .7 .7],'linewidth',1);
hold on;
plot3(fs2,real(fy),imag(fy),'color',[.5 .7 1],'linewidth',lw);
hold off;

axis([fs2(1) fs2(end) yLims2 yLims2]);
ax2.ColorOrder=min(origColorOrder+brighter,1);
ylabel('Re\{G\}');
zlabel('Imag');
xlabel('u');


subplot(324); % plot with colors

[scl,angImg]=imacx2(fAng');
    plot([fs2(i) fs2(i+1)],[fMag(i) fMag(i+1)],'color','k');
cla;
hold on;
for i=2:np-1
    color=max(0,min(1,angImg(i,:)));
%     color='w';
    plot([fs2(i) fs2(i+1)],[fMag(i) fMag(i+1)],'color',color,'linewidth',lw2);
end;
ysh=.014;
axis([fs2(1) fs2(end) -ysh 1+ysh]);
hold off;
xlabel('u');
ylabel('G(u)');
%%
ax3=subplot(326);
w=5;
line=zeros(np,np);
line(:,np/2-w:np/2+w)=repmat(fy,1,2*w+1);
pars.x=fs2;
pars.y=fs2;
imacx2(line,.5,pars);
ax3.YTick=[];
xlabel('u');


if writeImgs
    figName=['Gauss_shift.jpg'];
    disp(figName);
    print(figName,'-djpeg');
else
    pause(0.5);
end;



% function WriteOrPause(i,writeImgs,figNameExt)
% end

% % ax1=gca;
% set(gca,'colororder',min(origColorOrder+.3,1));
% % aVals=[1 2 3 5 8];
% % nt=numel(aVals);
% iTerms=zeros(numel(x2s),nt);
%
% cla;
% set(gca,'color',bkd);
% set(gca,'XColor','w');
% set(gca,'YColor','w');
% set(gca,'fontsize',18);
% ax1=gca;
% ax1.ColorOrderIndex=colorOrderIndex;
%
% for i=1:nt
%     a=aVals(i);
%     iTerms(:,i)=exp(-pi*(x2s/a).^2);
% end;

%     hold on;
%     baseColor=ax1.ColorOrder(ax1.ColorOrderIndex+i);
%     plot(x2s, terms(:,i),'-','color',termColor,'linewidth',1.5);
%     hold off
% end;

% %         for j=1:i
% %             text(x2s(1)+0.1,bases(1,j),num2str(fy(j+1),3),'verticalalignment','bottom',...
% %                 'color','w','fontsize',tfs);
% %         end;
%         hold off;
%         axis([-1 1 yLimits1]);
%         set(gca,'ytick',[]);
%         sy=sy+terms(:,i);
%
%
%         subplot(122); % ----------display sum of terms----------
%         hold on;
%         plot(x2s,y,'w--','linewidth',1.5,'markersize',8); % true function
%         %    plot(xs,[y sy],'linewidth',1.5);
%         %     if ~rectMode ||  mod(i,2)==1 % only if odd
%         plot(x2s,sy,'linewidth',lw);
%         %     end;
%         if i==nt
%             plot(x2s,sy(:,end),'color','w','linewidth',lw*1.33);
%         end;
%         hold off;
%         axis([-1 1 yLimits]);
%         WriteOrPause(i,writeImgs,figNameExt);
%     end;
%
% %%     Show the FT
%         subplot(122); % ------display individual terms---------
%         cla;
%         set(gca,'color',bkd);
%         set(gca,'XColor','w');
%         set(gca,'YColor',bkd);
%         set(gca,'fontsize',18);
%         wx=numel(xs)*16;
%         yx=Crop(y,numel(xs)*16); % pad by 16
%         fyx=real(fft(fftshift(yx)))/(w);
%         fyxs=Crop(fftshift(fyx),512);
%         xxs=(-256:255)'/128;
%         hold on;
%         plot(xxs,fyxs,'-','color',[1 .7 .7],'linewidth',2);
%         plot(xxs,xxs*0,'--','color',[.2 .7 .2]);
% %         plot(xs/2,y,'color',[1 .6 .6]);
%         for i=1:nt+1
%             plot((i-1)/8,2*fy(i),'wo','linewidth',lw*1.33);
%
%         end;
%         hold off;
%         axis([-2 2 1.2*min(fyxs) 1.2*max(fyxs)]);
%
%         WriteOrPause(0,writeImgs,['_ft_' figNameExt]);
% end;
%
