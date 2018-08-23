function m=f2CropRawMovie(m)
% Remove the strange pixels at the edge of a raw Falcon II movie
x1=97;
x2=4010;
y1=69;
y2=3989;

[nx,ny,nz]=size(m);
mask=zeros(nx,ny,'single');
mask(x1:x2,y1:y2)=1;
nmsk=sum(mask(:));
for i=1:nz
    m1=m(:,:,i);
    m1=mask.*m1;
    mean=sum(m1(:))/nmsk;
    m(:,:,i)=m1+(1-mask)*mean;
end;
