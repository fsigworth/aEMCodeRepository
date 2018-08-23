function [mxVals mxInds var cImage iptrs mxRimg totalVar]=rsVesicleCorrelation(m,mi,vesIndex,mbnOffset,partRadius,angleInds,eigenSet,noMask)
% function [mxVals mxInds var cimg]=rsVesicleCorrelation(m,mi,vesIndex,mbnOffset,partRadius,angleInds,eigenSet)
% Given an image m and mi structure, search for particles within the
% vesicle pointed to by vesIndex.  mbnOffset is the z-offset of the center
% of the membrane from the center of the original 3D model, in pixels with
% the same downsampling as m; partRadius (in
% m-pixels) is used for masking the cross-correlation function to the
% vicinity of the vesicle. angleInds is the 2d array of indices returned by
% rsListHemisphereAngles.
% The returned values are mxVals, the nt x nt 'image' of cross-correlation
% maxima; mxInds, the index giving gamma, hemi and io orientation of the
% particle; ptrs, the index of the projection direction (nt x nt x 2 array,
% where the third dimension is io); and cimg, the
% stack of eigenimage cross-correlations, ntxnt.
% The variance image (n x n) is returned as var.
%
% To determine the best-matching template at a given maximum of mxVals at
% image position (i,j), do this:
% g=mxInds(i,j);
% gamma=mod(g,nGamma);
% h=floor((g-1)/nGamma);
% hemi=mod(h,2)+1;  % 1 means upper hemisphere, 2 means lower
% ang=mxAngs(i,j);
% bestTemplate=xTemplates(:,:,ang,hemi,gamma);

if nargin<8
    noMask=0;
end;

n=size(m);

[nterms nHemi nGamma nAngs]=size(eigenSet.vList);
nio=2;
nRefs=nio*nHemi*nGamma;  % number of references to try at each pixel position
% nHemi is equal to 2.

maxDispersion=100;  % angstroms of radius
ds=mi.imageSize(1)/size(m,1);  % downsampling factor
ctr=[mi.vesicle.x(vesIndex) mi.vesicle.y(vesIndex)]/ds+1; % orig coords are 0-based.
r=mi.vesicle.r(vesIndex)/ds;
s0=mi.vesicle.s(vesIndex);
pixA=mi.pixA*ds;

% Determine nt, the size of the local cc function
ne=size(eigenSet.imgs,1);  % size of an image
nt=2*r+30/pixA+2*abs(mbnOffset)+ne... % maximum extent of a vesicle+particle
    +2*maxDispersion/pixA;  % plus surround
nt=2*ceil(nt*.55)+1;  % odd sized window for vesicle+particles, 10% oversized.
nt2=nt^2;


mt=ExtractImage(m,round(ctr),nt);  % Get the local part of the image
% ctrOffset=ctr-round(ctr);


% Get masks for the output cross-correlations
% msko=fuzzymask(nt,2,r+partRadius-mbnOffset,2);  % outside mask, extends
%                                     % partRadius beyond expected peak.
% mski=fuzzymask(nt,2,r+partRadius+mbnOffset,2);  % inside mask
if noMask
    msko=ones(nt,nt);
    mski=ones(nt,nt);
else
    extraR=2;
    msko=fuzzymask(nt,2,r-mbnOffset+extraR,2);  % outside mask, extends
    % partRadius beyond expected peak.
    mski=fuzzymask(nt,2,r+mbnOffset+extraR,2);  % inside mask
end;
% Get the larger mask for the variance calculation
if mbnOffset>0
    mskv=mski;
else
    mskv=msko;
end;
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
cimg=cimg.*repmat(msks,nHemi*nGamma,1)./(paArray);  % nRefs x npix
% cimg=cimg.*repmat(msks,nHemi*nGamma,1);  % nRefs x npix
cImage=zeros(nRefs,nt2);
cImage(:,pix)=cimg;
cImage=reshape(cImage,nRefs,nt,nt);

% find the maximum at each active pixel
[mxv mxi]=max(cimg);  % max over nRefs
% normalize by the chosen optimum match
for i=1:npix
    mxv(i)=mxv(i)/(s0*paArray(mxi(i),i));
end;
mximg=zeros(nt,nt);
mximg(pix)=mxv;      % image of maximum values, scaled

% Comput tIndexImg, giving the template index at each position in the image
evens=~mod(mxi,2);  % even reference number is inside-out
hemiAngs=pptrs(:,1);
hemiAngs(evens)=pptrs(evens,2);
tIndices=(ceil(mxi(:)/2)+max(0,hemiAngs-1)*nHemi*nGamma); % template indices
tIndexImg=zeros(nt,nt);
tIndexImg(pix)=tIndices;

lvars=squeeze(sum(cImage.^2)).*mskv;  % local variances computed from the ccs

mxVals=ExtractImage(mximg,round(ctr),n,1);
mxInds=ExtractImage(tIndexImg,round(ctr),n,1);
mxRimg=ExtractImage(rimg,round(ctr),n,1);
% mxAngs=ExtractImage(mxang,round(ctr),n,1);
var=ExtractImage(lvars,round(ctr),n,1);
totalVar=ExtractImage(squeeze(sum(ccs.^2,3)),round(ctr),n,1);