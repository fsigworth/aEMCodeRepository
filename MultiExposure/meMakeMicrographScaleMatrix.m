function [M,nOut,mi]=meMakeMicrographScaleMatrix(mi,ds)
% function [M,nOut,mi]=meMakeMicrographScaleMatrix(mi,ds)
% function M=meMakeMicrographScaleMatrix([],ds)
%  We create an affine transformation matrix M for the operation of
%  padding and downsampling a raw micrograph
%  to yield the output image size nOut, a 1x2 vector.
%  We make use of the mi fields imageSize and padImageSize
%  (If mi.padImageSize doesn't exist, it is initialized as
%     mi.padImageSize=NextNiceNumber(mi.imageSize,5,8)
%   so it is a multiple of 8 with largest prime factor 5. This is the only
%   reason for picking up the returned mi variable.)
%  The default downsampling factor ds is 1.
%  The returned matrix M is
%  [ds  0   -xsh
%   0   ds  -ysh
%   0   0    1 ]
%  Where the positive values [xsh ysh] are the shift of the origin in the
%  padded image, i.e. (mi.padImageSize-mi.imageSize)/2
%  The transformation from downsampled coordinates to the original coordinates will be
%  [x0; y0; 1]=M*[x; y; 1] or equivalently [x0 y0 1]=[x y 1]*M'
%  Note this the case of zero-based coordinates.
%  To map from original coordinates to downsampled coordinates,
%  [x;y;1]=inv(M)*[x0 y0 1];
%  or [x;y;1]=M\[x0 y0 1];
%  The corresponding crop-and-pad operation for the original micrograph is
%  mOut=Downsample(Crop(m0,mi.padSize),nOut);
%  --if the first argument is not a struct, we make a simple matrix M with
%  no shifts: [ds 0 0; 0 ds 0; 0 0 1].
% 
if numel(ds)<2 % isotropic downscaling
    ds=ds*[1 1];
end;
M=zeros(3,3,'single');
M(1,1)=ds(1);
M(2,2)=ds(2);
nOut=[];
if ~isa(mi,'struct') % No shifts required.
    return
else
    if ~isfield(mi,'padImageSize')
        mi.padImageSize=NextNiceNumber(mi.imageSize,5,8);
    end;
    M(:,3)=[round(mi.imageSize-mi.padImageSize)/2 1];
    nOut=round(mi.padImageSize./ds);
end;
