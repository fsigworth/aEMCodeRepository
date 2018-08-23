function [activeTrans,activeAlphasI]=reAssignActivePars(roi,thresholds)

% Handle translations.  Find translations above threshold, and expand by
% 1 pixel in the nt x nt array.
kernelT=true(3,3);  % kernel will expand by 1 pixel
pTrans=roi.pTrans;
nd=ndims(pTrans);
if nd<3  % see if we're handling pTrans as a vector
    [nt2,nim]=size(pTrans);
    nt=sqrt(nt2);
    pTrans=reshape(pTrans,nt,nt,nim);
end;

[nt,nt1,nim]=size(pTrans);
active=false(nt,nt,nim);
for i=1:nim
    a=pTrans(:,:,i)>thresholds(1);
    active(:,:,i)=BinaryConvolve(a,kernelT);
end;

if nd<3  % form active as a vector to match original pTrans
    active=reshape(active,nt2,nim);
end;

% Handle alphasI
nAI=size(roi.pAlphasI,1);
nA=nAI/r
a=roi.pAlphasI>thresholds(2);

