% IFTGaussian3.m
% Show ft pairs

writeImgs=1;
figNameExt='GaussFT.jpg';
index=2;
    
%     cd('/Users/fred/Documents/teaching/CMP710bCryoEM/710b-Videos/Video_lectures/2.3FT1D/IFT_figs/')
%     load origColorOrder % We'll modify the colors a little.
    
    bkd=[0 .14 .2]; % = 0 36 52 in Keynote.
    termColor=[.8 .8 .8];
    tfs=18;  % text label fontsize
    figure(1);
    clf;
    set(gcf,'color',bkd);
    set(gcf,'InvertHardCopy','off');

    %%     Show the FT
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
    plot(xxs,fy1xs,'-','color',[.5 .5 .5],'linewidth',2);
    plot(xxs,fyxs,'-','color',[1 .7 .7],'linewidth',2);
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
    
    
    WriteOrPause(index,writeImgs,figNameExt);
return

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