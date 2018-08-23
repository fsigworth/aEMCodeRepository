function [ampImg, mxInds, mxNCCImg]=spPickingCorrelation6(m,mi,eigenSet)
% function [mxVals mxInds mxNCCImg  mxValsV  mxRso  localVar]
%       =rspPickingCorrelation6(m,mi,vesIndex,mbnOffset,maskRadii,angleInds,eigenSet)
% Derived from rsVesicleCorrelation6.
% Given an image m and mi structure, search for particles within the
% vesicle pointed to by vesIndex.  mbnOffset is the z-offset of the center
% of the membrane from the center of the original 3D model, in pixels with
% the same downsampling as m; maskRadii are the radii of masks for masking
% the cross-correlation function to the vicinity of the vesicle, with radii
% for rso, iso and overall vesicle. angleInds is the 2d array of indices
% returned by rsListHemisphereAngles.
% Returned values:
% ampImg is the computed amplitude of the putative particle
%  at each pixel in the image.
% mxInds gives the index of the reference giving that value. This is the
%  index k into eigenSet.vList(:,k) or, alternatively,
%  vList(:,hemi,gamma,angs) where angs are enumerated according to
%  rsListHemisphereAngles.
% mxNCCImg is the maximum (over references) of the CC with normalized
%  references, such that ampImg=mxNCCImg./(ampList(mxInds)*s0), with s0
%  being the corresponding vesicle amplitude mi.vesicle.s(i,1).
%

[nterms, nHemi, nGamma, nAngs]=size(eigenSet.vList);
nRefs=nAngs*nHemi*nGamma;  % number of references to try at each pixel position

% Determine nt, the size of the local cc function
nt=size(m,1); %
nt2=nt^2;
paArray=reshape(eigenSet.ampList,nRefs,1);
%
% Compute the eigenimage cross-correlations
ccs=zeros(nt,nt,nterms);
fimg=fftn(m);
for i=1:nterms
    ccref=ifftshift(Crop(eigenSet.imgs(:,:,i),nt));  % pad each eigenimage
    ccs(:,:,i)=real(ifftn(fimg.*conj(fftn(ccref))));
    %     rimg=rimg+real(ifftn(fftn(ccs(:,:,i)).*fftn(ccref)))/nt;  % reconstructed image
end;
ccVec=reshape(ccs,nt2,nterms);

% Construct the stack of NCCs for each reference
nccPix=ccVec*reshape(eigenSet.vListNorm,nterms,nRefs);   % nt2 x nTerms * nTerms x nRefs

% find the maximum of the NCC at each active pixel in the masked CCs
[mxNCC, mxi]=max(nccPix,[],2);  % nPix x 1 : max over nRefs

% mxNCC is the maximum NCC value obtained over all the references.  The NCC
% is gotten by correlating normalized references with the original image.

% normalize by the chosen optimum match to get estimated amplitude
% If the image contains a*ref(j) at image position x,y then ampImg(x,y)
% will have the value a and mxInds(x,y) will be j.

ampPix=mxNCC./(paArray(mxi));

% Returned variables
mxInds=reshape(mxi,nt,nt);
mxNCCImg=reshape(mxNCC,nt,nt);
ampImg=reshape(ampPix,nt,nt);

