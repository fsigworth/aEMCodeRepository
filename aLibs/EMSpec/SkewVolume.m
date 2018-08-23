function v=SkewVolume(v0,zSlices,skew)

v=v0;
n=size(v,1);
nr=n/2;
for i=zSlices
    p=ToPolar(v0(:,:,i));
    for j=1:nr
        p(j,:)=circshift(p(j,:),[0,round(skew*j)]);
    end;
    v(:,:,i)=ToRect(p);
end;
