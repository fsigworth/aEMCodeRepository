function shift=AlignTransND(map,ref)


fmap=fftn(map-mean(map(:)));
fref=fftn(ref-mean(ref(:)));
cc=fftshift(real(ifftn(fmap.*conj(fref))));
n=size(cc);
ctr=fctr(n);
nd=ndims(map);
if nd==2 && any(n==1)
    nd=1;
end;
switch nd
    case 3
        [val,sh]=max3di(cc);
        shift=sh-ctr;
end;
