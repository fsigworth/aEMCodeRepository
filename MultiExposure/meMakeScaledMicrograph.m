function [mOut,M]=meMakeScaledMicrograph(m0,mi,ds)
% Pad and normalize a raw micrograph, downsampling by ds if desired.
% The transform matrix M is returned. (See meMakeMicrographScaleMatrix for
% details.)
if numel(ds)==1
    ds=[ds ds];
end;
mNorm=mi.imageNormScale*(m0-mi.imageMedian);
[M,nOut]=meMakeMicrographScaleMatrix(mi,ds);
mOut=Downsample(Crop(mNorm,mi.imageSize),nOut);
% (We could have applied the normalization after downsampling...)
