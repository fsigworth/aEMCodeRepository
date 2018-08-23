function fVols=reFourierInsertion(projs,degAngles)
% function fVols=reFourierInsertion(projs,degAngles)
% projs should have size n x n x nAngs x 2 x nSets
% where projs(:,:,:,1) are the class means and projs(:,:,:,2) are the
% norms.
sz=size(projs);
if numel(sz)<5
    sz(5)=1;
end;
n=sz(1);
nRefs=sz(3);
nVols=prod(sz(4:end));  
projs=reshape(projs,n,n,nRefs,nVols);

fVols=gridMakeNullFT(n,3);
fVols(nVols,1)=fVols;

angles=rsDegToEuler(degAngles);


parfor iVol=1:nVols
% for iVol=1:nVols
    p=projs(:,:,:,iVol);
    fvol=gridMakeNullFT(n,3);
    
    for i=1:nRefs
        nslice=gridMakePaddedFT(p(:,:,i));
        fvol=gridInsertPlane(nslice,fvol,angles(i,:));
    end
    fVols(iVol)=fvol;
end;
fVols=reshape(fVols,sz(4:end));
