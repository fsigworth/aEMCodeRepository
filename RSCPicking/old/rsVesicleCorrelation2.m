function [mxVals mxInds var cimg]=rsVesicleCorrelation2(m,mi,vesIndex,mbnOffset,partRadius,angleInds,eigenSet)
% function [mxVals mxInds cimg]=rsVesicleCorrelation(m,mi,vesIndex,mbnOffset,partRadius,angleInds,eigenSet)
% Given an image m and mi structure, search for particles within the
% vesicle pointed to by vesIndex.  mbnOffset is the z-offset of the center
% of the membrane from the center of the original 3D model; partRadius (in
% pixels) is used for masking the cross-correlation function to the
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

n=size(m);
mskBorder=2;

[nterms nAngs nHemi nGamma]=size(eigenSet.vList);
nio=2;
% nHemi is equal to 2.

maxDispersion=100;  % angstroms of radius
ds=mi.imageSize(1)/size(m,1);  % downsampling factor
ctr=[mi.vesicle.x(vesIndex) mi.vesicle.y(vesIndex)]/ds+1; % orig coords are 0-based.
r=mi.vesicle.r(vesIndex)/ds;
s0=mi.vesicle.s(vesIndex);

% Determine nt, the size of the local cc function
ne=size(eigenSet.imgs,1);  % size of an image
nt=2*r+30/mi.pixA+2*abs(mbnOffset)+ne... % maximum extent of a vesicle+particle
    +2*maxDispersion/mi.pixA;  % plus surround
nt=NextNiceNumber(nt*1.1,7,1);  % largest prime factor is 7, may be odd.
% nt=2*ceil(nt*.55)+1;  % odd sized window for vesicle+particles, 10% oversized.
nt2=nt^2;
mt=ExtractImage(m,round(ctr),nt);  % Get the local part of the image
ctrOffset=ctr-round(ctr);

% Compute the cross-correlations
ccs=zeros(nt,nt,nterms);
fimg=fftn(mt);

for i=1:nterms
    ccref=ifftshift(Crop(eigenSet.imgs(:,:,i),nt));  % pad each eigenimage
    ccs(:,:,i)=real(ifftn(fimg.*conj(fftn(ccref))));
end;

rad=Radius(nt);

msko=rad<(r+partRadius-mbnOffset+mskBorder);  % outside mask, extends
                                    % partRadius beyond expected peak.
pixo=find(msko);
mski=rad<r+partRadius+mbnOffset+mskBorder;
pixi=find(mski);

% 
ptrsR=rsGetAngleMap(nt,r-mbnOffset,0,angleInds);  % outside-out map
ptrsI=rsGetAngleMap(nt,r+mbnOffset,1,angleInds);  % inside-out map
ptrsRso=ptrsR(pixo);
ptrsIso=ptrsI(pixi);

% These reshapings take ~300 ms
% This following loop takes 54% of the total time in this function!!!
% vArray is nGamma*4 x nterms x nt^2
%  gamma changes most rapidly, then upper/lower, then rso/iso
for i=1:2  % rso/iso
        vArrayRso=eigenSet.vList(:,ptrsRso,:,:);  % (nterms pixels nHemi nGamma)
        vArrayIso=eigenSet.vList(:,ptrsIso,:,:);  % iso vesion
    end;
end;

% projamps=reshape(eigenSet.ampList,nAngs,2*nGamma);
paArray=zeros(nt^2,4*nGamma);
for i=1:2
    for j=1:2
        k=(j-1)+2*(i-1);
        paArray(:,nGamma*k+1:nGamma*(k+1))...
            =squeeze(eigenSet.ampList(ptrs(:,:,i),j,:));
    end;
end;

% a little faster:
% paArray=[squeeze(eigenSet.ampList(ptrso,1,:)),...
%          squeeze(eigenSet.ampList(ptrso,2,:)),...
%          squeeze(eigenSet.ampList(ptrsi,1,:)),...
%          squeeze(eigenSet.ampList(ptrsi,2,:))];  % nt^2 x 4nGamma
% 
ccArray=shiftdim(reshape(ccs,nt*nt,nterms),1);   % now nterms x nt^2

% Construct the stack of cross-correlations in factor space
cimg=zeros(4*nGamma,nt*nt);
vArray(nterms,npix,etc.) * ccArray(npix,nterms)
    cimg(:,i)=vArray(:,:,i)*ccArray(
for i=1:nt*nt
    cimg(:,i)=vArray(:,:,i)*ccArray(:,i);
end;

% Mask the cross-correlations outside the vesicle
msko=fuzzymask(nt,2,r+partRadius-mbnOffset,2);  % outside mask, extends
                                    % partRadius beyond expected peak.
mski=fuzzymask(nt,2,r+partRadius+mbnOffset,2);  % inside mask
cimg=cimg.*[repmat(msko(:)',nGamma*2,1); repmat(mski(:)' ,nGamma*2,1)]./(s0*paArray');
% Get the mask for the variance calculation
if mbnOffset>0
    mskv=mski;
else
    mskv=msko;
end;

% Normalize so the nominal amplitude is 1,and maximize over gamma and io
[mxv mxi]=max(cimg);  % W'll make 2D images of max values
pts=reshape(ptrs,nt*nt,2);  % get the pointers as 2 vectors
mxa=pts(:,1);
mxa(mxi>(2*nGamma))=pts(mxi>(2*nGamma),2); % select the 2nd vector for io

mxang=reshape(mxa,nt,nt);
mximg=reshape(mxv,nt,nt);
mxind=reshape(mxi,nt,nt);
lvars=sum(ccs.^2,3).*mskv;  % local variances computed from the ccs
cimg=reshape(cimg,nGamma*4,nt,nt);  % Direct cross-correlation

mxVals=ExtractImage(mximg,round(ctr),n,1);
mxInds=ExtractImage(mxind,round(ctr),n,1);
mxAngs=ExtractImage(mxang,round(ctr),n,1);
var=ExtractImage(lvars,round(ctr),n,1);
