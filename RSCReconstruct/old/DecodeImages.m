% DecodeImages
nta=size(transImgs,2);
qta=zeros(32,32,nta);
for i=1:nta
    q=zeros(32,32);
    q(pix)=transImgs(:,i);
    qta(:,:,i)=q;
end;
qr=zeros(32,32,nRefs);
for i=1:nRefs
    q=zeros(32,32);
    q(pix)=ctRefs(:,i,iMi);
    qr(:,:,i)=q;
end;
