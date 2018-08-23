function rings=makeRings(imgSize,nRings)
[x y]=meshgrid([1:imgSize],[1:imgSize]);
center=imgSize/2+0.5;
rI=sqrt((x-center).^2+(y-center).^2);
rings=zeros([imgSize imgSize nRings]);
dr=imgSize/(2*(nRings));
for i=1:nRings,
    tmp= (rI>=(i-1)*dr) & (rI<i*dr);
    rings(:,:,i)=tmp;
end