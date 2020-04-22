% IFTGaussianToDelta.m
% Make a movie of limiting Gaussians
% The movie is written out as a series of jpeg images.
% fs Nov 2018
writeImgs=1;
% nameEnd='_square.jpg';
figNameExt='_gauss.jpg';


bkd=[0 .14 .2]; % = 0 36 52 in Keynote.
termColor=[.8 .8 .8];
brighter=0.3;
lw=2;
tfs=18;  % text label fontsize
fs=18;
figure(1);
clf;
set(gcf,'color',bkd);
set(gcf,'InvertHardCopy','off');
np=200;
w=np/4;
xs=(-2*w:2*w-1)'/w;
x2s=xs/2;


subplot(121); % ------display individual terms---------
cla;
% set(gca,'colororder',min(origColorOrder+brighter,1));
set(gca,'color',bkd);
set(gca,'XColor','w');
set(gca,'YColor','w');
set(gca,'fontsize',18);
ax1=gca;
ax1.ColorOrderIndex=colorOrderIndex;

subplot(122); % for the FTs
cla;
% set(gca,'colororder',min(origColorOrder+brighter,1));
set(gca,'color',bkd);
set(gca,'XColor','w');
set(gca,'YColor','w');
set(gca,'fontsize',18);
ax2=gca;
ax2.ColorOrderIndex=colorOrderIndex;


aVals=[1 2 4 8 16];
nt=numel(aVals);
terms=zeros(numel(x2s),nt);
iTerms=zeros(numel(x2s),nt);
yLims1=[0 4];
yLims2=[0 1.3];

for i=1:nt
    a=aVals(i);
    terms(:,i)=a*exp(-pi*(a*x2s).^2);
    iTerms(:,i)=exp(-pi*(x2s/a).^2);
end;

xt=.8;
for i=1:nt
    subplot(121);
    plot(x2s, terms(:,1:i),'linewidth',lw);
    axis([-1 1 yLims1]);
    ax1.ColorOrder=min(origColorOrder+brighter,1);

    subplot(122); % plot the FTs
    plot(x2s, iTerms(:,1:i),'linewidth',lw);
    hold on;
    for j=1:i
        if j==nt
%             xt=0.5;
            str=['a=' num2str(aVals(j)) '   ' ];
            alString='right';
        else
            str=[num2str(aVals(j))];
            alString='left';
        end;
    text(xt,exp(-pi*(xt/aVals(j))^2),str,'color','w', ...
        'fontsize', fs,'verticalalignment','bottom','horizontalalignment',alString);
    end;
    hold off;
    axis([-1 1 yLims2]);
    ax2.ColorOrder=min(origColorOrder+brighter,1);

drawnow;
if writeImgs
    figName=['delta_' num2str(i) figNameExt];
    disp(figName);
    print(figName,'-djpeg');
else
    pause(0.5);
end;
    
    
    
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
