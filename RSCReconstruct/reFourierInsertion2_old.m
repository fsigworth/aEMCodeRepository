function fVols=reFourierInsertion2(projs,refAngles)
% function fVols=reFourierInsertion(projs,ri)
% projs should have size n x n x nAngs x nNorm x nVols
% We use ri to specify projection angles, and ri.symmetry to do multiple
% insertions.
% where projs(:,:,:,1,:) are the class means and projs(:,:,:,2,:) are the
% norms.  The returned fVols is an array of structures, 2 x nVol in size.

sz=size(projs);
if numel(sz)<5  % number of volumes to be determined overall.
    nVols=1;
end;
n=sz(1);
nAngs=sz(3);
nNorm=sz(4);
nV2=nNorm*nVols;

projs=reshape(projs,n,n,nAngs,nV2);

fVols1=gridMakeNullFT(n,3);
for i=1:nV2
    fVols(i,1)=fVols2;  % array of structures
end;

% Get the reference angles in the Euler system, radians
% [betas, gammas]=reGetRefAngles(ri,refInds);
refAnglesE=rsDegToEuler(refAngles);

if nAngs ~= size(refAngles,1)
    error(['Wrong number of ref angles', num2str(nAngs)]);
end;

if mod(ri.symmetry,2)==0
    sym2=ri.symmetry/2;
elseif ri.symmetry==1
    sym2=1;
else
    error(['Unsupported symmetry: ' num2str(ri.symmetry)]);
end;

% Get the extra gamma rotation for symmetric insertions
refAngleOffsetE=rsDegToEuler([0 0 180/sym2]);

parfor iVol=1:nV2
% for iVol=1:nVols
    p=projs(:,:,:,iVol);
    fvol=gridMakeNullFT(n,3);
    
    for i=1:nAngs
        nslice=gridMakePaddedFT(p(:,:,i));
        for j=1:sym2
            angles=refAnglesE(i,:)+(j-1)*refAngleOffsetE;
            fvol=gridInsertPlane(nslice,fvol,angles);
        end;
    end
    fVols(iVol)=fvol;
end;
fVols=reshape(fVols,sz(4:end));
