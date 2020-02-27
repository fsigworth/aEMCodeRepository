function SubplotLabel(x,y,label,fontsize,color)
% function SubplotLabel(x,y,label,fontsize,color)
% Given normalized coordinates, put a label in the current axes.
% The default fontsize is 18, and default color is [0 0 0] (=black).
% These values are persistent, so you don't have to give them after the
% first call.
% If label is numeric, then label=1 -> A, label=2 -> B etc.
% Color can be either a char, e.g. 'r', a scalar (gives gray level, 0 to 1),
% or a three element vector (RGB), where [1 1 1] is white.
%
% (Note that you must not use inf in setting any of the limits of the axes.)

persistent fsz clr

if numel(fsz)<1 % first time through, initialize
    fsz=18;
    clr=[0 0 0];
end;

if isnumeric(label)
    label=char(64+label);
end;

if nargin<4
    fontsize=fsz;
end;
if nargin<5
    color=clr;
end;

if isnumeric(color) && numel(color)<3
    color=color(1)*[1 1 1];
end;

h=gca;
xLim=h.XLim;
yLim=h.YLim;
u=xLim(1)+x*(xLim(2)-xLim(1));
v=yLim(1)+y*(yLim(2)-yLim(1));
text(u,v,label,'fontsize',fontsize,'color',color);
