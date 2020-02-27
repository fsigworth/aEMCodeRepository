function vlSetDarkGraphics(fontSize,bkdColor,doReset)
% function vlSetDarkGraphics(fontSize,bkdColor,doReset)
% Set the dark graphics standard for video lectures.
% doReset restores to the normal values.
% Defaults are fontSize=14, bkdColor=dark blue, doReset=0.

if nargin<3
    doReset=0;
end;
if nargin<2 || numel(bkdColor)<1
    bkdColor=[0 .14 .2]; % = 0 36 52 in Keynote.
end;
if nargin<1 || numel(fontSize)<1
    fontSize=18;  % default value.
end;

if ~doReset
    bkdColor=[0 .14 .2]; % = 0 36 52 in Keynote.
    set(groot,'defaultFigureColor',bkdColor);
    set(groot,'defaultFigureInvertHardCopy','off');
    set(groot,'defaultAxesColor',bkdColor);
    set(groot,'defaultAxesXColor','w');
    set(groot,'defaultAxesYColor','w');
    set(groot,'defaultAxesFontSize',fontSize);
    set(groot,'defaultLineLineWidth',2);
    
else
    rmv='remove';
    set(groot,'defaultFigureColor',rmv);
    set(groot,'defaultFigureInvertHardCopy',rmv);
    set(groot,'defaultAxesColor',rmv);
    set(groot,'defaultAxesXColor',rmv);
    set(groot,'defaultAxesYColor',rmv);
    set(groot,'defaultAxesFontSize',rmv);
    set(groot,'defaultLineLineWidth',rmv);
end;
