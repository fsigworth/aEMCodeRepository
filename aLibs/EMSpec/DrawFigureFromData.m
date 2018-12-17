function DrawFigureFromData(dat,nr,nc,ipanel)
% function DrawFigureFromData(dat,nr,nc,ipanel)
% Given a figure data structure dat, draw a single figure panel or a
% multi-panel figure; the latter starts with subplot index ipanel.
% Only the first argument is required.
% 
% Supported are only images and simple line plots. For a single panel
% dat is a struct such as
%     dat.xs=[-10 10];  % optional
%     dat.ys=[-10 10];  % optional, but note that xs and ys must precede the image field
%     dat.image=randn(20,20);
%     dat.xLabel='X-label';
%     dat.yLabel='Y-label';
%     dat.title='Title';
% Or, alternatively,
%     dat.xplot=1:20;
%     dat.yplot=randn(20,4);
% with the same label options as above.
% to create an entire figure, use nested fields e.g.
%     mdat.thePlot.yplot=randn(10,2);
%     mdat.thePlot.title='Random plot';
%     mdat.theImage.image=randn(10,10);
%     mdat.theImage.title='Random image';
%     DrawFigureFromData(mdat);  % draws two panels in subplot(1,2,.) figure.

if nargin<4
    ipanel=1;
end;
if nargin<3
    nc=0;
end;
if nargin<2
    nr=0;
end;

graphicsNames={'image' 'ploty'};
names=fieldnames(dat);
singleAxes=0;
for i=1:numel(graphicsNames)
    if any(strcmp(graphicsNames{i},names))
        singleAxes=1;
        break;
    end;
end;

if singleAxes
    nr=max(nr,1);
    nc=max(nc,1);
    subplot(nr,nc,ipanel);
    DrawPanel(dat);
else
    nPanels=numel(names)+ipanel-1;  % total subpanels needed
    if nr*nc<nPanels  % not enough subplots
        nr=floor(sqrt(nPanels));
        nc=ceil(nPanels/nr);
    end;
    for i=1:numel(names)
        subplot(nr,nc,ipanel+i-1)
        if isa(dat.(names{i}),'struct')
            DrawPanel(dat.(names{i}));
        end;
    end;
end;
end

function DrawPanel(s)
        xs=[];
        ys=[];
        imagePresent=0;
        plotx=[];
    fn=fieldnames(s);
    for i=1:numel(fn);
        switch fn{i}
            case 'exec'
                eval(s.eval);
            case 'xs'
                xs=s.xs;
                if numel(ys)<1
                    ys=xs;
                end;
            case 'ys'
                ys=s.ys;
            case 'image'
                imags(xs,ys,s.image);
                imagePresent=1;
            case 'plotx'
                plotx=s.plotx;
            case 'ploty'
                if imagePresent
                    hold on;
                end;
                if numel(plotx)>0
                    plot(plotx,s.ploty);
                else
                    plot(s.ploty);
                end;
                if imagePresent
                    hold off;
                end;
            case 'xLabel'
                xlabel(s.xLabel);
            case 'yLabel'
                ylabel(s.yLabel);
            case 'title'
                title(s.title,'interpreter','none');
            case 'titleTex'
                title(s.titleTex);
            case 'axis'
                 axis(s.axis)
        end;

    end;
end