function [ampImg, mxInds, mxNCCImg, mxValsV, mxRso, localVar]=rsVesicleCorrelation5(m,mi,vesIndex,mbnOffset,maskRadii,angleInds,eigenSet)
% function [mxVals mxInds mxNCCImg  mxValsV  mxRso  localVar]
%       =rsVesicleCorrelation(m,mi,vesIndex,mbnOffset,maskRadii,angleInds,eigenSet)
% Given an image m and mi structure, search for particles within the
% vesicle pointed to by vesIndex.  mbnOffset is the z-offset of the center
% of the membrane from the center of the original 3D model, in pixels with
% the same downsampling as m; maskRadii are the radii of masks for masking
% the cross-correlation function to the vicinity of the vesicle, with radii
% for rso, iso and overall vesicle. angleInds is the 2d array of indices
% returned by rsListHemisphereAngles.
% Returned values: mxVals is the maximum NCC value at each pixel, and
% mxInds gives the index of the reference giving that value.  filtImg is
% the image reconstructed from the eigenimages. mxValsV is the same as
% mxVals except it isn't masked.  mxRso tells whether the best-matching
% particle orientation is RSO (=1) or ISO (=0).  localVar is the variance
% computed from the first few terms of the eigenimage expansion.

% To determine the best-matching template at a given maximum of mxVals at
% image position (i,j), do this:
% g=mxInds(i,j);
% gamma=mod(g,nGamma);
% h=floor((g-1)/nGamma);
% hemi=mod(h,2)+1;  % 1 means upper hemisphere, 2 means lower
% ang=mxAngs(i,j);
% bestTemplate=xTemplates(:,:,ang,hemi,gamma);

[nterms, nHemi, nGamma, nAngs]=size(eigenSet.vList);
nVarTerms=min(1,nterms);  % number of terms to sum for local variance.
nio=2;
nRefs=nio*nHemi*nGamma;  % number of references to try at each pixel position
% nHemi is equal to 2.

ds=mi.imageSize(1)/size(m,1);  % downsampling factor
ctr=[mi.vesicle.x(vesIndex) mi.vesicle.y(vesIndex)]/ds+1; % orig coords are 0-based.
r=mi.vesicle.r(vesIndex)/ds;
s0=mi.vesicle.s(vesIndex);
% % pixA=mi.pixA*ds;

% Determine nt, the size of the local cc function
% ne=size(eigenSet.imgs,1);  % size of an image
nt=NextNiceNumber(ceil(2*max(maskRadii))+1,7,1); % largest factor is 7
nt2=nt^2;


mt=ExtractImage(m,round(ctr),nt);  % Get the local part of the image
% ctrOffset=ctr-round(ctr);


% Get masks for the output cross-correlations
% msko=fuzzymask(nt,2,r+partRadius-mbnOffset,2);  % outside mask, extends
%                                     % partRadius beyond expected peak.
% mski=fuzzymask(nt,2,r+partRadius+mbnOffset,2);  % inside mask
    msko=fuzzymask(nt,2,maskRadii(1),.5)>.5;  % outside-out mask
    % partRadius beyond expected peak.
    mski=fuzzymask(nt,2,maskRadii(2),.5)>.5;  % inside-out mask
    mskv=fuzzymask(nt,2,maskRadii(3),.5)>.5;  % overall vesicle mask (largest)

% Use mskv to find the active pixels
pix=find(mskv(:)>.5);  % column vector
npix=numel(pix);

% Look up the angle indices for each active pixel
iptrs=zeros(nt,nt,2);
iptrs(:,:,1)=rsGetAngleMap(nt,r-mbnOffset,0,angleInds);  % outside-out map
iptrs(:,:,2)=rsGetAngleMap(nt,r+mbnOffset,1,angleInds);  % inside-out map
ptrs=reshape(iptrs,nt^2,2);
pptrs=ptrs(pix,:); % npix x 2, pointers for each active pixel

vArrayNorm=zeros(nterms,nio,nHemi,nGamma,npix);
for i=1:2  % rso/iso
%     vArray(:,i,:,:,:)=eigenSet.vList(:,:,:,pptrs(:,i));
    vArrayNorm(:,i,:,:,:)=eigenSet.vListNorm(:,:,:,pptrs(:,i));
end;
% vArray=reshape(vArray,nterms,nRefs,npix);
vArrayNorm=reshape(vArrayNorm,nterms,nRefs,npix);

paArray=zeros(nio,nHemi,nGamma,npix);
for i=1:2
    paArray(i,:,:,:)=eigenSet.ampList(:,:,pptrs(:,i));
end;
paArray=reshape(paArray,nRefs,npix);

% Compute the eigenimage cross-correlations
ccs=zeros(nt,nt,nterms);
% rimg=zeros(nt,nt);
fimg=fftn(mt);
for i=1:nterms
    ccref=ifftshift(Crop(eigenSet.imgs(:,:,i),nt));  % pad each eigenimage
    ccs(:,:,i)=real(ifftn(fimg.*conj(fftn(ccref))));
%     rimg=rimg+real(ifftn(fftn(ccs(:,:,i)).*fftn(ccref)))/nt;  % reconstructed image
end;
ccVec=reshape(ccs,nt2,nterms);
ccPix=ccVec(pix,:);  % npix x nterms
localVarPix=sum(ccPix(:,1:nVarTerms).^2,2);  % npix x 1: local variance at each pixel
ccPixNorm=ccPix;  % ** no division by sigma_img

% Construct the stack of NCCs for each reference
nccPix=zeros(nRefs,npix);
for i=1:npix
    nccPix(:,i)=(ccPixNorm(i,:)*vArrayNorm(:,:,i))';
end;

% Special for only outside-out particles....
% Apply the masks
nccPixMasked=min(nccPix(:))*ones(nRefs,npix,'single');
pto=find(msko(pix));
pti=find(mski(pix));
for i=1:2:nRefs
    nccPixMasked(i,pto)=nccPix(i,pto);
    nccPixMasked(i+1,pti)=nccPix(i+1,pti);
end;


% nccPix=reshape(nccPix,nio,nHemi*nGamma,npix);
% nccPixMasked=-inf*ones(nio,nHemi*nGamma,npix,'single');
% 
% nccPixMasked(1,:,msko(pix))=nccPix(1,:,msko(pix)); % nRefs x npix
% % nccPixMasked(2,:,mski(pix))=nccPix(2,:,mski(pix)); %.*shiftdim(repmat(mski(pix)',nHemi*nGamma,1),-1);
% % % Don't allow iso points to be chosen
% % mskiPix=repmat(~mski(pix)',nHemi*nGamma,1);
% % nccPixMasked(2:2:nRefs,:)=nccPixMasked(2:2:nRefs,:)-1e9*mskiPix; % nearly infinite
% % nccPixMasked=nccPix;

% find the maximum of the NCC at each active pixel in the masked CCs
% % nccPixMasked=reshape(nccPixMasked,nRefs,npix);
[mxNCC, mxi]=max(nccPixMasked);  % 1 x nPix : max over nRefs

% normalize by the chosen optimum match to get estimated amplitude
denormNCC=mxNCC';
ampPix=zeros(npix,1);
for i=1:npix
    ampPix(i)=denormNCC(i)/(s0*paArray(mxi(i),i));
end;
ampImg=zeros(nt,nt);
ampImg(pix)=ampPix;      % image of maximum values, scaled

% Comput tIndexImg, giving the template index at each position in the image
evens=~mod(mxi,2);  % even reference number (from masked cc) is inside-out
hemiAngs=pptrs(:,1);
hemiAngs(evens)=pptrs(evens,2);
tIndices=(ceil(mxi(:)/2)+max(0,hemiAngs-1)*nHemi*nGamma); % template indices
tIndexImg=zeros(nt,nt);
tIndexImg(pix)=tIndices;

% lvars=squeeze(sum(cImage.^2)).*mskv;  % local variances computed from the ccs

% Returned variables
% % mxVals=ampImg;
% ampImgV=zeros(nt,nt);
% ampImgV(pix)=ampVPix;      % image of maximum values, scaled
mxInds=tIndexImg;
mxNCCImg=zeros(nt,nt);
mxNCCImg(pix)=mxNCC;
localVar=zeros(nt,nt);
localVar(pix)=localVarPix;

% mxValsV=ampImgV;% *** removed
mxValsV=ampImg; 
mxRso=uint8(zeros(nt,nt));  % image where 1=right side out, 0= ISO.
mxRso(pix)=~evens;
