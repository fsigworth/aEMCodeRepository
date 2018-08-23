function mi=rsCreateVesicleModel(mi,thk,rise,overshoot)
if nargin<4
    overshoot=0;
end;

% membrane model
vLipid=1.6;
% Create the model, which is sampled in units of the original pixel size.
nm0=ceil((thk/2+2*rise)/mi.pixA)*2+1;  % array for vesicle model
mi.vesicleModel=fuzzymask(nm0,1,thk/mi.pixA/2,rise/mi.pixA)...
    *vLipid;  % units of V.A per voxel
for j=-1:2:1
    p=round((nm0+1+j*(thk-rise/2)/mi.pixA)/2);
mi.vesicleModel(p)=mi.vesicleModel(p)*(1+overshoot);
end;
