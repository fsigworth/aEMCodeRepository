function [fves]=meMakeVesicleReferences(n,mbnModel,pixA,rPars)
rmin=rPars(1)/pixA;
rmax=rPars(2)/pixA;
rstep=rPars(3)/pixA;
nrsteps=ceil((rmax-rmin)/rstep)+1;

fves=single(complex(zeros(n,n,nrsteps)));
for i=1:nrsteps
    r=rmin+(i-1)*rstep;  % radius in pixels
    v = VesicleFromModel(n,r,mbnModel);
    fves(:,:,i)=fftn(ifftshift(v));    % FT of vesicle at origin.
end;
