function SubplotLabel(x,y,label,fontsize,color)
% function SubplotLabel(x,y,label,fontsize,color)
% Given normalized coordinates, put a label in the current axes.
% The default fontsize is 18, and default color is [0 0 0] (=black).

if nargin<4
    fontsize=18;
end;
if nargin<5
    color=[0 0 0];
end;
if numel(color)<3
    color=color(1)*[1 1 1];
end;

h=gca;
xLim=h.XLim;
yLim=h.YLim;
u=xLim(1)+x*(xLim(2)-xLim(1));
v=yLim(1)+y*(yLim(2)-yLim(1));
text(u,v,label,'fontsize',fontsize,'color',color);
