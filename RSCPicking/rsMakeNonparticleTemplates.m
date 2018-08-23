function templates=rsMakeNonparticleTemplates(n,nGamma,mi)
% Make a set of 2xnGamma non-particle templates.  There are nGamma rotated
% edges of vesicles, and nGamma spherical blobs.

% Find a typical vesicle
vesRadius=median(mi.vesicle.r);
vesAmp=median(mi.vesicle.s);
[err vindex]=min(abs(mi.vesicle.r-vesRadius)/vesRadius...
                +abs(mi.vesicle.s-vesAmp)/vesAmp);
            
% multiply the vesicle model by the voxel size.
vesModel=meDownsampleVesicleModel(mi.vesicleModel,ds)*ds*mi.pixA;
r0=round(vesRadius/vs);
n=NextNiceNumber(4*vesRadius/ds);
v=VesicleFromModel(n,r0,vesModel);
for i=1:nGamma
    alpha=(i-1)*360/nGamma;
    ctr=n/2+1;
    x0=sind(alpha)*r0+ctr;
    y0=cosd(alpha)*r0+ctr;
    

