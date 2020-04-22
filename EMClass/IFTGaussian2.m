% IFTGaussian4.m
% Version 4 has computing the FT integral

writeImgs=0;
% nameEnd='_square.jpg';
nameEnd='_short2.jpg';

for rectMode=1
    
    
%     cd('/Users/fred/Documents/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.3FT1D/IFT_figs/')
%     load origColorOrder % We'll modify the colors a little.
    
    bkd=[0 .14 .2]; % = 0 36 52 in Keynote.
    termColor=[.8 .8 .8];
    tfs=18;  % text label fontsize
    figure(1);
    clf;
    set(gcf,'color',bkd);
    set(gcf,'InvertHardCopy','off');
    np=200;
    w=np/4;
    xs=(-2*w:2*w-1)'/w;
    x2s=xs;
    
    subplot(122);
    cla;
    set(gca,'color',bkd);
    set(gca,'XColor','w');
    set(gca,'YColor','w');
    
    set(gca,'fontsize',18);
    ax2=gca;
    
    if rectMode
        nt=10;
        y=abs(xs)<.5;
%         y(51)=.5;
%         y(151)=.5;
        termThresh=.02;
        yscl=.5;
        yLimits=[-.4 1.4];
        yLimits1=[1.3  nt/2+.5];
        lw=1;
        set(gca,'colororder',min(1,origColorOrder+.1));
        colorOrderIndex=5;
        figNameExt=['_rect' nameEnd];
        marker='w--';
        fcnText='y=rect(x)';
        tPos=[.3 1.1];
        efs=tfs; % equation fontsize
    else
        nt=6;
        y=exp(-pi*xs.^2);
        termThresh=0;
        yscl=1;
        yLimits=[-.2 1.1];
        yLimits1=[.8 nt-2];
        lw=1.5;
        set(gca,'colororder',min(origColorOrder+.1,1));
        figNameExt=['_gauss' nameEnd];
        colorOrderIndex=9;
        marker='w--';
        fcnText='y=e^{-\pi x^2}';
        fcnText2='g(x) = e^{-\pi x^2}';
        tPos=[.6 .8];
        efs=tfs+8;
    end;
    fy=real(fft(fftshift(y)))/(2*w);
    fy(abs(fy)<termThresh)=0;
    
    % plot the original function
    hold on;
    plot(x2s,y,'w-','linewidth',1.5,'markersize',8); % true function
    sy=fy(1)/2+0*xs;
    % plot(x2s,sy,'linewidth',1);
    hold off;
    axis([-2 2 yLimits]);
    text(tPos(1),tPos(2),fcnText,'color','w','fontsize',efs);
    
    
    WriteOrPause(0,writeImgs,figNameExt);
    cla
    set(gca,'colororder',min(origColorOrder+.3,1));
    
    subplot(121);
    ax1=gca;
    terms=zeros(np,nt);
    base=zeros(1,nt);
    for i=1:nt
        subplot(121); % ------display individual terms---------
        cla;
        set(gca,'color',bkd);
        set(gca,'XColor','w');
        set(gca,'YColor','w');
        set(gca,'fontsize',18);
        ax1=gca;
        ax1.ColorOrderIndex=colorOrderIndex;
        terms(:,i)=cos(i/2*pi*xs)*fy(i+1);
        if rectMode
            base(i)=nt/2-.35*i+.25*(i==1);
        else
            base(i)=nt/2-.35*i+2*abs(fy(i+1));
        end;
        bases=repmat(base(1,1:i),np,1);
        hold on;
        %     baseColor=ax1.ColorOrder(ax1.ColorOrderIndex+i);
        plot(x2s,bases,'--','color',[.8 .8 .4]);
        plot(x2s, yscl*terms(:,1:i)+bases,'-','color',termColor,'linewidth',1.5);
        for j=1:i
            text(x2s(1)+0.1,bases(1,j),num2str(2*fy(j+1),3),'verticalalignment','bottom',...
                'color','w','fontsize',tfs);
        end;
        hold off;
        axis([-2 2 yLimits1]);
        set(gca,'ytick',[]);
        sy=sy+terms(:,i);
        
        
        subplot(122); % ----------display sum of terms----------
        hold on;
        plot(x2s,y,'w--','linewidth',1.5,'markersize',8); % true function
        %    plot(xs,[y sy],'linewidth',1.5);
        %     if ~rectMode ||  mod(i,2)==1 % only if odd
        plot(x2s,sy,'linewidth',lw);
        %     end;
        if i==nt
            plot(x2s,sy(:,end),'color','w','linewidth',lw*1.33);
        end;
        hold off;
        axis([-2 2 yLimits]);
        WriteOrPause(i,writeImgs,figNameExt);
    end;
    return
%%     Show the FT
clf;
set(gcf,'color',bkd);
xRange=4;
uRange=2;
showPoints=0;
wr=50;
xs=(-xRange*wr:xRange*wr-1)'/wr;
          y=(exp(-pi*(xs-2).^2)+exp(-pi*(xs+2).^2));
         y1=2*exp(-pi*xs.^2);
subplot(121);
    cla;
    set(gca,'color',bkd);
    set(gca,'XColor','w');
    set(gca,'YColor','w');
    
    set(gca,'fontsize',18);
    fy=real(fft(fftshift(y)))/wr;
    fy(abs(fy)<termThresh)=0;
    
    % plot the original function
    hold on;
    plot(xs,y,'w-','linewidth',1.5,'markersize',8); % true function
    sy=fy(1)/2+0*xs;
    % plot(x2s,sy,'linewidth',1);
    hold off;
    axis([-xRange xRange 0 1.19]);
%     text(tPos(1),tPos(2),fcnText2,'color','w','fontsize',efs);
   ylabel('g(x)');
    xlabel('x');

        subplot(122); % ------display the FT---------
        cla;
        set(gca,'color',bkd);
        set(gca,'XColor','w');
        set(gca,'YColor','w');
        set(gca,'fontsize',18);
        wf=50; % expansion of FT
        df=wf*2*xRange;
        xxs=(-uRange*df:uRange*df-1)'/(df);
        wx=numel(xs)*wf;

        yx=Crop(y,numel(xs)*wf); % pad by wf
        fyx=real(fft(fftshift(yx)))/wf;
        fyxs=Crop(fftshift(fyx),numel(xxs));
        
        y1x=Crop(y1,numel(xs)*wf); % pad by wf
        fy1x=real(fft(fftshift(y1x)))/wf;
        fy1xs=Crop(fftshift(fy1x),numel(xxs));
        
        hold on;
        plot(xxs,fyxs,'-','color',[1 .7 .7],'linewidth',2);
        plot(xxs,fy1xs,'-','color',[.5 .5 .5],'linewidth',2);
        plot(xxs,xxs*0,'--','color',[.2 .7 .2]);
%         plot(xs/2,y,'color',[1 .6 .6]);
if showPoints
    for i=1:nt+1
            plot((i-1)/(2*uRange),fy(i),'wo','markersize',10,'linewidth',2);
        end;
end;
hold off;
        axis([-uRange uRange 1.2*min(fyxs) 1.2*max(fyxs)]);
           ylabel('G(u)');
    xlabel('u');

    
        WriteOrPause(0,writeImgs,['_ft_' figNameExt]);
end;

function WriteOrPause(i,writeImgs,figNameExt)
drawnow;
if writeImgs
    figName=[num2str(i) figNameExt];
    disp(figName);
    print(figName,'-djpeg');
else
    pause(0.5);
end;
end