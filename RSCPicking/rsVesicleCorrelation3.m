function [mxVals mxInds filtImg mxValsV mxRso]=rsVesicleCorrelation3(m,mi,vesIndex,mbnOffset,maskRadii,angleInds,eigenSet)
% function [mxVals mxInds filtImg mxValsV mxRso]
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
% particle orientation is RSO (=1) or ISO (=0).

% To determine the best-matching template at a given maximum of mxVals at
% image position (i,j), do this:
% g=mxInds(i,j);
% gamma=mod(g,nGamma);
% h=floor((g-1)/nGamma);
% hemi=mod(h,2)+1;  % 1 means upper hemisphere, 2 means lower
% ang=mxAngs(i,j);
% bestTemplate=xTemplates(:,:,ang,hemi,gamma);

[nterms nHemi nGamma nAngs]=size(eigenSet.vList);
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
    msko=fuzzymask(nt,2,maskRadii(1),.5);  % outside-out mask
    % partRadius beyond expected peak.
    mski=fuzzymask(nt,2,maskRadii(2),.5);  % inside-out mask
    mskv=fuzzymask(nt,2,maskRadii(3),.5);  % overall vesicle mask (largest)

% Use mskv to find the active pixels
pix=find(mskv(:)>.5);  % column vector
npix=numel(pix);

% Look up the angle indices for each active pixel
iptrs=zeros(nt,nt,2);
iptrs(:,:,1)=rsGetAngleMap(nt,r-mbnOffset,0,angleInds);  % outside-out map
iptrs(:,:,2)=rsGetAngleMap(nt,r+mbnOffset,1,angleInds);  % inside-out map
ptrs=reshape(iptrs,nt^2,2);
pptrs=ptrs(pix,:); % npix x 2, pointers for each active pixel

vArray=zeros(nterms,nio,nHemi,nGamma,npix);
for i=1:2  % rso/iso
    vArray(:,i,:,:,:)=eigenSet.vList(:,:,:,pptrs(:,i));
end;
vArray=reshape(vArray,nterms,nRefs,npix);

paArray=zeros(nio,nHemi,nGamma,npix);
for i=1:2
    paArray(i,:,:,:)=eigenSet.ampList(:,:,pptrs(:,i));
end;
paArray=reshape(paArray,nRefs,npix);

% Compute the eigenimage cross-correlations
ccs=zeros(nt,nt,nterms);
rimg=zeros(nt,nt);
fimg=fftn(mt);
for i=1:nterms
    ccref=ifftshift(Crop(eigenSet.imgs(:,:,i),nt));  % pad each eigenimage
    ccs(:,:,i)=real(ifftn(fimg.*conj(fftn(ccref))));
    rimg=rimg+real(ifftn(fftn(ccs(:,:,i)).*fftn(ccref)))/nt;  % reconstructed image
end;

ccVec=reshape(ccs,nt2,nterms);
ccPix=ccVec(pix,:);  % npix x nterms

% Construct the stack of cross-correlations for each reference
cimg=zeros(nRefs,npix);
for i=1:npix
    cimg(:,i)=(ccPix(i,:)*vArray(:,:,i))';
end;
% mask the cross-correlations and convert back to image
msks=[msko(pix) mski(pix)]';  % 2 x npix
cimgM=cimg.*repmat(msks,nHemi*nGamma,1)./(paArray);  % nRefs x npix
msksV=mskv(pix)';
cimgV=cimg.*repmat(msksV,nRefs,1)./paArray;  % whole-vesicle cross-correlation

% find the maximum at each active pixel
[mxv mxi]=max(cimgM);  % max over nRefs

% normalize by the chosen optimum match
% mxvUnNorm=mxv.*localSD;
mxvUnNorm=mxv;  %%%%%%%%%%%%%%%%%?????

for i=1:npix
    mxv(i)=mxvUnNorm(i)/(s0*paArray(mxi(i),i));
end;
mximg=zeros(nt,nt);
mximg(pix)=mxv;      % image of maximum values, scaled

% Do the same for the whole-vesicle image.
[mxV mxiV]=max(cimgV);
for i=1:npix
    mxV(i)=mxV(i)/(s0*paArray(mxiV(i),i));
end;
mximgV=zeros(nt,nt);
mximgV(pix)=mxV;      % image of maximum values, scaled

% Comput tIndexImg, giving the template index at each position in the image
evens=~mod(mxi,2);  % even reference number is inside-out
hemiAngs=pptrs(:,1);
hemiAngs(evens)=pptrs(evens,2);
tIndices=(ceil(mxi(:)/2)+max(0,hemiAngs-1)*nHemi*nGamma); % template indices
tIndexImg=zeros(nt,nt);
tIndexImg(pix)=tIndices;

% lvars=squeeze(sum(cImage.^2)).*mskv;  % local variances computed from the ccs

% Returned variables
mxVals=mximg;
mxInds=tIndexImg;
filtImg=rimg;
mxValsV=mximgV;

mxRso=uint8(zeros(nt,nt));  % image where 1=right side out, 0= ISO.
mxRso(pix)=~evens;
