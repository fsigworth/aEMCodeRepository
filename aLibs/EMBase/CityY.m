function pvals=CityY(vals)
if length(vals)==numel(vals) % a 1d array
    vals=vals(:);
    pvals=repmat(vals',2,1);
    pvals=pvals(:);
else
    vals=squeeze(vals); % better be 2d
    sz=size(vals);
    pvals=repmat(shiftdim(vals,-1),2,1,1);
    pvals=reshape(pvals,sz(1)*2,sz(2));
end;

