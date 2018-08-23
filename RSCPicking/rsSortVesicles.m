function miNew=rsSortVesicles(mi)
% Sort the vesicle locations as a raster starting from the top left.
% The fields vesicle.x, y, r, s and optionally ok are all sorted.

nBands=16;  % number of horizontal bands

xs=mi.vesicle.x;
ys=mi.vesicle.y;

y1=mi.imageSize(2);
dy=y1/nBands;  % make that many horizontal bands
allInds=[];    % We'll put the sorted indices here.
for y0=y1-dy:-dy:0
    p=find(ys>=y0+1 & ys<y0+dy);
    if numel(p)>0
        bandxs=mi.vesicle.x(p);
        [sxs indx]=sort(bandxs);
        allInds=[allInds p(indx)];
    end;
end;
miNew=mi;
miNew.vesicle.x=mi.vesicle.x(allInds);
miNew.vesicle.y=mi.vesicle.y(allInds);
miNew.vesicle.r=mi.vesicle.r(allInds);
miNew.vesicle.s=mi.vesicle.s(allInds);
if isfield(miNew.vesicle,'ok')
    miNew.vesicle.ok=mi.vesicle.ok(allInds);
end;
