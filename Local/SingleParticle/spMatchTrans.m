function [aliImgs dxy]=spMatchTrans(stack,sp,refs)
% function aliImgs=spMatchTrans(stack,sp,refs)
% Perform a single-reference translational alignment of the stack images to
% the already-assigned class reference.  Don't update the sp structure.

[n ny nim]=size(stack);
dxy=zeros(nim,2);
aliImgs=single(zeros(n,ny,nim));

for i=1:nim
    j=sp.class(i);
if j>0
    im=stack(:,:,i);
    fim=fftn(im);
    fref=fftn(refs(:,:,j));
    cc=fftshift(real(ifftn(fim.*conj(fref))));
    [mx x y]=max2di(cc);
    dxy(i,:)=[x-n/2-1 y-ny/2-1];
    P=FourierShift([n ny],-dxy(i,:));
    im2=real(ifftn(P.*fim));
    aliImgs(:,:,i)=im2;
end;
end;

