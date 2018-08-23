function pvals=CityX(xvals,yvals)
% Modify x-coordinates for making a cityscape plot.
% Usage:
% xs=
if nargin<2
    yvals=xvals;
end;
if size(yvals,1)==1
    yvals=shiftdim(yvals,1);
end;
xvals=xvals(:);
pvals=repmat(xvals',2,1);
pvals=pvals(:);
pvals(1)=[];
if size(yvals,1)==numel(xvals)-1
    pvals(end)=[];
else
    pvals(end+1)=pvals(end);
end;

end

