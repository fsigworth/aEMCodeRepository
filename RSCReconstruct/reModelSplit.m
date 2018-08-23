function mo1=reModelSplit(moi,indices)
% Copy fields of moi into the reduced structre mo1.  indices gives the
% indices of imgAmps and pRefs to copy; it can also be a vector of booleans.
mo1.sigmaN=moi.sigmaN;
mo1.sigmaC=moi.sigmaC;
mo1.sigmaG=moi.sigmaG;
mo1.b0=moi.b0;
mo1.pVols=moi.pVols;
mo1.a=moi.a;
mo1.imgAmps=moi.imgAmps(indices);
mo1.activeTrans=moi.activeTrans(:,:,indices);
if numel(moi.pRefs)>1
    mo1.pRefs=moi.pRefs(indices,:);
else
    mo1.pRefs=moi.pRefs;
end;
