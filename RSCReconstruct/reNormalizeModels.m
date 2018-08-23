function vModel=reNormalizeModels(ri,vRecon)
% function vModel=reNormalizeModels(ri,vRecon)
% Rescale the volume(s) to match the original one.

nVols=size(vRecon,4);
vModel=vRecon;
for iVol=1:nVols
    v=vRecon(:,:,:,iVol);
    vsd=sqrt(v(:)'*v(:)/numel(v));  % Get the global SD of all vols.
    vModel(:,:,:,iVol)=v*ri.volSD(iVol)*ri.refScale/vsd;
end;
