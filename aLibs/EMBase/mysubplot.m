function mysubplot(nr,nc,ind,xs,ys,xo,yo)
% function mysubplot(nr,nc,ind,xs,ys,xo,yo)
% Like subplot() but leaves narrow gaps between axes.
if nargin<2  % handle 3 digits like '221'
    nr0=nr;
    nr=floor(nr0/100);
    nr0=nr0-100*nr;
    nc=floor(nr0/10);
    nr0=nr0-10*nc;
    ind=nr0;
end;
if nargin<5
    xs=.03;
    ys=.025;
end;
if nargin<6
    xo=0;
    yo=0;
end;
xEdgeL=.05+xo;
xEdgeR=.03;
yEdgeL=.05+yo;
yEdgeU=.03;
yEdgeT=0;
hpos=mod(ind-1,nc);  % zero-based horizontal index
vpos=nr-1-floor((ind-1)/nc); % zero-based vertical index
xwid=(1-xEdgeL-xEdgeR-(nc-1)*xs)/nc;
ywid=(1-yEdgeL-yEdgeU-(nr-1)*ys-yEdgeT)/nr;
xorg=xEdgeL+xs*hpos+hpos*xwid;
yorg=yEdgeL+ys*vpos+vpos*ywid;

pos=[xorg yorg xwid ywid];
subplot('position',pos);
% Turn off interior tick labels
if hpos>0
    set(gca,'YTickLabel','');
end;
if vpos>0
    set(gca,'XTickLabel','');
end;
% imags(randn(128))
