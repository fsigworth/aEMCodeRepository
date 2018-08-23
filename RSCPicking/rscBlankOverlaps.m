function [ccMax, ovMask]=rscBlankOverlaps(mi,ccMax,overlapRadius)
% Blank regions of the cross-correlation map where two vesicles overlap,
% or come within overlapRadius (given in angstroms).
% if overlapRadius==0, nothing is masked.
% 
n=size(ccMax);
ds=mi.imageSize(1)/n(1);
% Convert everything to pixels on ccMax
vesx=mi.vesicle.x/ds+1;
vesy=mi.vesicle.y/ds+1;
vesr=mi.vesicle.r/ds;
ovr=overlapRadius/(mi.pixA*ds);
numv=numel(mi.vesicle.x);

ovMask=single(zeros(n));  % We'll mask the CC where msk > 1 due to overlap.
if overlapRadius==0
    return
end;

for i=1:numv
    if all(mi.vesicle.ok(i,1:2:3)) % vesicle exists and was fitted
    r1=vesr(i)+ovr;  % r1 = radius + ovr of the ith vesicle
%     Clip the vesicle coordinates
    x0=max(1,floor(vesx(i)-r1));
    x1=min(n(1),ceil(vesx(i)+r1));
    y0=max(1,floor(vesy(i)-r1));
    y1=min(n(2),ceil(vesy(i)+r1));
    [lx,ly]=ndgrid(x0:x1,y0:y1);
%     Mask a disc of radius r1
    ovMask(x0:x1,y0:y1)=ovMask(x0:x1,y0:y1)+...
        ((lx-vesx(i)).^2+(ly-vesy(i)).^2 <= r1^2);
    end;
end;
% mask the cc function
ccMax=ccMax.*(ovMask<1.5);
