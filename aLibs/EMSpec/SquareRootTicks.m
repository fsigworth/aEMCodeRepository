function SquareRootTicks(xvals,num,yaxis)
% function SquareRootTicks(xvals,reciprocal)
% Put x-axis ticks in the appropriate place for plotting square-root
% values.
% e.g. plot(sqrt(freqs),spect); SquareRootTicks(freqs);
% if nargin<2
% reciprocal=0;
% end;
if nargin<2
    num=20;
end;
if nargin<3
    yaxis=false;
end;
xmin=min(xvals);
xmax=max(xvals);
dx=Step125(xmax/num);
dx2=Step125(dx/4);
dx4=Step125(dx2/4);
x0=floor(xmin/dx)*dx;
tickVals=[x0:dx4:dx2-dx4 dx2:dx2:dx-dx2 dx:dx:xmax+dx-eps]';
% tickVals
% if reciprocal
%     vals=1/tickVals;
% else
    vals=tickVals;
% end;
labels=num2str(vals);
if yaxis
    set(gca,'ytick',sqrt(tickVals),'yticklabel',labels);
else
    set(gca,'xtick',sqrt(tickVals),'xticklabel',labels);
end;
