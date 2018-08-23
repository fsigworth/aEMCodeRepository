function mRefs=reMakeMatchedRefs(refs,si,ri,roi)
% From the information in the roi structure, compute the best-matching
% reference image for each member of the image stack, appropriately
% transformed, scaled and ctf-filtered. Subtracting this from the
% original should give pure noise.

nim=numel(si.miIndex);
n=size(si.ctfs,1);

sOffset=(sqrt(size(roi.pTrans,1))+1)/2;
mRefs=zeros(n,n,nim,'single');

% for i=1:nim
for i=1:nim
    [~,iRef,iVol]=max2d(roi.pRefs(:,:,i));
    T=roi.imgTA(:,i);
    if all(T>0)
        tRef=circshift(mean(refs(:,:,iRef,iVol,:),5),T(1:2)'-sOffset);
        rRef=rsRotateImage(tRef,ri.alphasI(T(3)));
        cRef=real(ifftn(fftn(rRef).*ifftshift(si.ctfs(:,:,si.miIndex(i)))));
        mRefs(:,:,i)=roi.imgAmps(i)*cRef;
    end;
end;

