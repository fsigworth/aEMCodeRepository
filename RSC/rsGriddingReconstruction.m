function [vol norm volf]=rsGriddingReconstruction(angles,imgs,norms,k)
% Perform 3D reconstruction using a stack of projection images and a
% corresponding stack of norms, both real-space images.  The (nangs x 3) list of
% angles is in degrees; k is the Wiener constant.  For Frealign-like reconstruction
% the images should be pre-filtered by the ctf, and the norms should look
% approximately like delta functions in the center of the images.  Hence
% one image each in the stack of imgs and norms might be obtained as
%
%  c=ifftshift(CTF(n,pixA,ctPars));  % shift frequency origin to edge
% imgs(:,:,i)=real(ifftn(fftn(img).*c));% raw image filtered 'again' by ctf
% norms(:,:,i)=fftshift(real(ifftn(c.^2)));
% The size of the images must be a multiple of 4.
%
doParallel=1;

angs=rsDegToEuler(angles);  % convert to Euler angles in radians.

ks=3;
[n ny nangs]=size(imgs);

comp=gridMakePreComp(n,ks);

normVol=gridMakeNullFT(n,3);  % Volume to receive normalization
bpVol=normVol;                % back-projected volume

% snorms=shiftdim(norms,2);

if doParallel
    
    % We do Fourier insertions into two volumes.
     nworkers=4;
%     tic
    nlangs=ceil(nangs/nworkers);
    nvsum=normVol.PadFT;
    bpsum=bpVol.PadFT;
    parfor j=1:nworkers
        i0=(j-1)*nlangs+1;
        i1=min(nangs, j*nlangs);
        nv=gridMakeNullFT(n,3);  % Volume to receive normalization
        bv=gridMakeNullFT(n,3);  % Volume to receive normalization
        for i=i0:i1
            fslice0=gridMakePaddedFT(norms(:,:,i));
            nv=gridInsertPlane(fslice0,nv,angs(i,:));
            fslice=gridMakePaddedFT(imgs(:,:,i));
            bv=gridInsertPlane(fslice,bv,angs(i,:));
        end;
        
        nvsum=nvsum+nv.PadFT;
        bpsum=bpsum+bv.PadFT;
    end;
    
    normVol.PadFT=nvsum;
    bpVol.PadFT=bpsum;
    
    
%     disp('par done');
%     toc
    
else % non-parallel
    for i=1:nangs
        fslice0=gridMakePaddedFT(norms(:,:,i));
        normVol=gridInsertPlane(fslice0,normVol,angs(i,:));
    end;
    
    for i=1:nangs
        fslice=gridMakePaddedFT(imgs(:,:,i));
        bpVol=gridInsertPlane(fslice,bpVol,angs(i,:));
    end;
end;

fmsk=fuzzymask(bpVol.np1,3,bpVol.np*0.45,n*0.05);

% Wiener filter "normalization"
volf=bpVol;
volf.PadFT=bpVol.PadFT.*fmsk./(k^2+normVol.PadFT);

% Get the reconstructed volume.
vol=gridRecoverRealImage(volf,comp);
if nargout>1
    norm=gridRecoverRealImage(normVol,comp);
end;
