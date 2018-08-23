function mi=meInsertVesicleModel(mi,thickA,riseA)
% function mi=meInsertVesicleModel(mi,thickA,riseA)
% Insert a generic membrane model into the mi structure. The optional
% arguments thickA are riseA are the thickness and the border width, in
% angstroms.

vLipid=1.6;  % inner potential in volts

if nargin<2
    thickA=50;
end;
if nargin<3
    riseA=6;
end;

% Create the model, which is sampled in units of the original pixel size.
nm0=ceil(30/mi.pixA)*2+1;  % array for vesicle model; 60A nominal
mi.vesicleModel=fuzzymask(nm0,1,thickA/mi.pixA/2,riseA/mi.pixA)...
    *vLipid*mi.pixA;  % units of V.A per voxel
mi.vesicleModelSampInt=mi.pixA;

