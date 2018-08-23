function map=deHotPixelMap(img,sds)
% function map=deHotPixelMap(img,sds)
% Create a map where ones indicate outlier pixels.  An outlier is a pixel 
% where the Laplacian is greater than sds * standard deviation of the Laplacian
% of the entire image.  
% To correct an image typically you would afterwards perform
%  RemoveOutliers(img,map);

[nx ny]=size(img);
locMeans=single(zeros(nx,ny));

for i=-1:2:1
    for j=-1:2:1
        locMeans=locMeans+circshift(img,[i j]);
    end;
end;
locMeans=locMeans/4;
d=img-locMeans;
sd=std(d(:));
map=abs(d)>sds*sd;

