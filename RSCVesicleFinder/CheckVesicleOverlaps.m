% CheckVesicleOverlaps
ves=mi.vesicle;
nv=numel(mi.vesicle.x);
n=mi.imageSize/8;
vSum=zeros(n,'single');
vBad=zeros(n,'single');
overlap=zeros(nv,1);
for i=1:nv
    ok=all(mi.vesicle.ok(i,1));
    v0=meMakeModelVesicles(mi,n,i,0,0);
    v1=double(v0<.8*min(v0(:)));
    overlap(i)=(v1(:)'*vSum(:))/(v1(:)'*v1(:));
if overlap(i)<.3 && ok
    vSum=vSum+v1;
elseif ok
    vBad=vBad+v1;
end;
end;

return
%%
n=960;
for i=1:nv
    v0=meMakeModelVesicles(mi,n,i,0,0);
    msk=v0<.8*min(v0(:));
    imags(msk)
    drawnow;
end;