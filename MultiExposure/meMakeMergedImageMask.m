function msk=meMakeMergedImageMask(n,Tmats, edgeWidth, weights)
% function msk=meMakeMergedImageMask(n,Tmats, edgeWidth, weights)
% Create a mask of size n corresponding to images merged according to the
% affine matrices Tmats. Masked region of msk is zeros.

[t1 t2 nim]=size(Tmats);
if nargin<3 || all(edgeWidth==0)
    edgeWidth=n/512;  % same criterion as in meCombineImages.
end;
edgeWidth=ceil(edgeWidth);
if nargin<4
    weights=ones(1,nim);
end;

n2x=NextNiceNumber(1.5*n(1)); % padded mask

mfree=ones(n-2*edgeWidth);  % free area of a masked image
msk=Crop(mfree,n);      % Add masked borders
mx=Crop(msk,n2x);        % padded version

for i=2:nim
    if weights(i)
        T=Tmats(:,:,i);
        T(1:2,3)=T(1:2,3).*(n./n2x)';  % scale down the shifts according to the padding
        mxT=AffineTransform(mx,T,'nearest');
        msk=msk.*Crop(mxT,n);
    end;
end;
