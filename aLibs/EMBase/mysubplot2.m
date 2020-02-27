function mysubplot2(nr,nc,ind,ll,sp,ru)
% function mysubplot(nr,nc,ind,ll,sp,ru)
% Like subplot() but allows narrow gaps between axes.
% The extra arguments are 1x2 vectors giving [x y] coordinates:
% ll 'left, lower' is the origin of the bottom left axes;
% sp 'spacing' is the spacing in x and y
% ru 'right, upper is the upper right corner
%
% if you want to suppress internal tick labels, after drawing say
% set(gca,'XTickLabel','');

persistent org space fin labels initialized

if nargin<2  % handle 3 digits like '221'
    nr0=nr;
    nr=floor(nr0/100);
    nr0=nr0-100*nr;
    nc=floor(nr0/10);
    nr0=nr0-10*nc;
    ind=nr0;
end;
if numel(initialized)==0 % put in default values
    disp('initializing mysubplot');
    org=[.06 .06];
    space=[.03 .03];
    fin=[.02 .02];
    initialized=1;
end;

if nargin<4
    ll=org;
else
    org=ll;
end;
if nargin<5
    sp=space;
else
    space=sp;
end;
if nargin<6
    ru=fin;
else
    fin=ru;
end;
if nargin<7
    intLabels=labels;
elseif numel(intLabels)<2
    intLabels=[1 1]*intLabels;
    labels=intLabels;
else
    labels=intLabels;
end;

hpos=mod(ind-1,nc);  % zero-based horizontal index
vpos=nr-1-floor((ind-1)/nc); % zero-based vertical index

span=(1-ll-ru+sp)./[nc nr]-sp;
edge=ll+[hpos vpos].*(sp+span);

pos=[edge span];

subplot('position',pos); % we assume normalized coordinates
