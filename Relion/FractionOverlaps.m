function [m,frac]=FractionOverlaps(mi);
    thk=4;
    minRadius=100;
    maxRadius=400;
    ds=8;  % downsampling from original image
    nv=numel(mi.vesicle.x);
    n=mi.padImageSize/ds;
    
     mi.vesicle.ok(:,2)=mi.vesicle.ok(:,2) & ...
         mi.vesicle.s(:,1)>0.5*median(mi.vesicle.s(:,1)) & ...
        mi.vesicle.r(:,1)>minRadius & mi.vesicle.r(:,1)<maxRadius;
    goodVes=(mi.vesicle.ok(:,2));
    goodVesInds=find(goodVes);
    if numel(goodVesInds)>0
        m=meMakeModelVesicleDiscs(mi,n,goodVesInds,thk);
        m=m.*meGetMask(mi,n);
        frac=sum(m(:)>1)/sum(m(:)>0.5);
    else
        m=zeros(n,'single');
        frac=0;
    end;
end
%         imags(v1);
%         drawnow;
        
