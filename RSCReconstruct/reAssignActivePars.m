function [activeTrans,activeAlphas,activeRVs]=reAssignActivePars(ri,roi)
% Generalization of reAssignActiveTrans to handle active Alphas as well.
% Checks for p(Trans) or P(AlphasI) smaller than ri.thresholds(1,2),
% respectively.

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
activeTrans=false(nt,nt,nim);
for i=1:nim
    a=pTrans(:,:,i)>ri.thresholds(1);
    activeTrans(:,:,i)=BinaryConvolve(a,kernelT);
end;

if nd<3  % form active as a vector to match original pTrans
    activeTrans=reshape(activeTrans,nt2,nim);
end;

% Handle alphasI
nA=size(ri.alphas,1);
nI=size(ri.alphas,2);
a=roi.pAlphas>ri.thresholds(2);  % above threshold
% Expand by one at each end
%a=reshape(a,nA,nI,nim);
a(2:nA,:,:)=a(2:nA,:,:) | a(1:end-1,:,:);
a(1:nA-1,:,:)=a(1:nA-1,:,:) | a(2:nA,:,:);
activeAlphas=a;

% Handle refs
a=roi.pRefs>ri.thresholds(3);
activeRVs=BinaryConvolve(a,true(3,3));


