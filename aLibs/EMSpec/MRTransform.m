function alignedimages=MRTransform(imgs,xytfr,ds,msk)
% function alignedimages=MRTransform(imgs,xytfr,ds)
% Given the transform variables xytfr from MRAlign and the optional downsampling factor
% ds (>1 if xytfr was obtained on downsampled images), shift and rotate
% imgs to give alignedimages.
if nargin<3
    ds=1;
end;
[n, ny, nim]=size(imgs);
if nargin<4
    msk=fuzzymask(n,2,n*0.45,n*.2);
end;
alignedimages=single(zeros(n,n,nim));
for j=1:nim
    theta=-xytfr(j,3);
    if theta==0
        alignedimages(:,:,j)=circshift(imgs(:,:,j),-round(xytfr(j,1:2)*ds));
    else
        alignedimages(:,:,j)=grotate(msk.*circshift(imgs(:,:,j),...
            -round(xytfr(j,1:2)*ds)),-xytfr(j,3));
    end;
end;
