function fVols=reFourierInsertion(projs,refAngles,symmetry,useParFor,doDisp)
% function fVols=reFourierInsertion(projs,ri,symmetry)
% projs are n-D and typically have size n x n x nAngs x nNorm x nVols
% where nNorm = 1 or 2 to denote classes or norm projections.
% If symmetry>2 we do insertions with multiple gamma values.
% The returned fVols is an n-D array of structures with the dimensions
% being dimensions nNorm x nVols x (any other dimensions of projs)

if nargin<3
    symmetry=1;
end;
if nargin<4
    useParFor=0;
end;
refTime=0;
if nargin<5
    doDisp=0;
else
    tic;
end;
sz=size(projs);
if numel(sz)<5
    sz(end+1:5)=1; % make sure elements 4... are a matrix
end;
n=sz(1);
nAngs=sz(3);
nVtotal=prod(sz(4:end));  % total Fourier volumes to create

projs=reshape(projs,n,n,nAngs,nVtotal);

fVols1=gridMakeNullFT(n,3);
fVols=fVols1;
for i=2:nVtotal
    fVols(i,1)=fVols1;  % array of structures
end;

% Get the reference angles in the Euler system, radians
% [betas, gammas]=reGetRefAngles(ri,refInds);
refAnglesE=rsDegToEuler(refAngles);

if nAngs ~= size(refAngles,1)
    error(['Wrong number of ref angles', num2str(nAngs)]);
end;

doVolSym=false;
if mod(symmetry,2)==0  % symmetry is even
    if symmetry<=4
        sym2=1;
        doVolSym=true;
    else
        sym2=symmetry/2;
    end;
elseif symmetry==1
    sym2=1;
else
    error(['Unsupported symmetry: ' num2str(ri.symmetry)]);
end;

% Get the extra gamma rotation for symmetric insertions
refAngleOffsetE=rsDegToEuler([0 0 360/symmetry])-rsDegToEuler([0 0 0]);

if useParFor
    parfor iVol=1:nVtotal
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
        if doVolSym
            fvol.PadFT=Symmetrize3(fvol.PadFT,symmetry);
        end;
        fVols(iVol)=fvol;
    end;
else
    if doDisp 
        disp(['volume | angle of [' num2str([nVtotal nAngs]) ']' ]);
    end;
    for iVol=1:nVtotal
        p=projs(:,:,:,iVol);
        fvol=gridMakeNullFT(n,3);
        
        for i=1:nAngs
        if doDisp && toc > doDisp
            tic;
            disp([iVol i]);
        end;
            nslice=gridMakePaddedFT(p(:,:,i));
            for j=1:sym2
                angles=refAnglesE(i,:)+(j-1)*refAngleOffsetE;
                fvol=gridInsertPlane(nslice,fvol,angles);
            end;
        end
        if doVolSym
            fvol.PadFT=Symmetrize3(fvol.PadFT,symmetry);
        end;
        fVols(iVol)=fvol;
    end;
end;

fVols=reshape(fVols,sz(4:end));
