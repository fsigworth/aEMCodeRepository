function fVols=reFourierInsertion2(projs,refAngles,symmetry,useParFor)
% function fVols=reFourierInsertion(projs,ri,symmetry)
% projs are n-D and typically have size n x n x nAngs x nNorm x nVols
% where e.g. nNorm = 1 or 2 to denote classes or norm projections.
% If symmetry>1 we do insertions with multiple gamma values.
% The returned fVols is an n-D array of structures with the dimensions
% nNorm x nVols x (any other dimensions of projs)

if nargin<3
    symmetry=1;
end;
if nargin<4
    useParFor=0;
end;

sz=size(projs);
if numel(sz)<5
    sz(end+1:5)=1; % make sure elements 4... are a matrix
end;
n=sz(1);
nAngs=sz(3);
nVols=prod(sz(4:end));  % total Fourier volumes to create

projs=reshape(projs,n,n,nAngs,nVols);

fVols1=gridMakeNullFT(n,3);
fVols=fVols1;
for i=2:nVols
    fVols(i,1)=fVols1;  % array of structures
end;

% Get the reference angles in the Euler system, radians
% [betas, gammas]=reGetRefAngles(ri,refInds);
refAnglesE=rsDegToEuler(refAngles);

if nAngs ~= size(refAngles,1)
    error(['Wrong number of ref angles', num2str(nAngs)]);
end;

% Get the extra gamma rotation for symmetric insertions
refAngleOffsetE=rsDegToEuler([0 0 360/symmetry])-rsDegToEuler([0 0 0]);

if useParFor
    jMax=nVols*symmetry;
    tVols=fVols1;  % temporary partial volumes
    for j=2:jMax
        tVols(j,1)=fVols1;
    end;
    parfor j=1:jMax
        iVol=ceil(j/symmetry);
        iSym=mod(j,symmetry);
        p=projs(:,:,:,iVol);
        for i=1:nAngs
            nslice=gridMakePaddedFT(p(:,:,i));
            angles=refAnglesE(i,:)+(iSym-1)*refAngleOffsetE;
            tVols(j)=gridInsertPlane(nslice,tVols(j),angles);
        end;
    end;
    for j=1:jMax
        iVol=ceil(j/symmetry);
        fVols(iVol).PadFT=fVols(iVol).PadFT+tVols(j).PadFT;
    end;
else        
    for iVol=1:nVols
        p=projs(:,:,:,iVol);
        
        for i=1:nAngs
            nslice=gridMakePaddedFT(p(:,:,i));
            for j=1:symmetry
                angles=refAnglesE(i,:)+(j-1)*refAngleOffsetE;
                fVols(iVol)=gridInsertPlane(nslice,fVols(iVol),angles);
            end;
        end
    end;
end;

fVols=reshape(fVols,sz(4:end));
