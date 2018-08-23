function aliImgs=TransformImages(imgs, xytfr)
% function aliImgs=TransformImages(imgs, xytfr);
% Transform the image stack imgs according to the shifts, rotates and flips
% in the nim x 5 array xytfr.
% xytfr contains the coordinate system operations needed to return the
% image to align to the reference.

nim=size(imgs,3);
nim2=size(xytfr,1);

if nim2<nim
    error('Not enough elements in xytfr');
end;

theta=xytfr(:,3);
iflip=xytfr(:,4);
im=shiftf(imgs,-xytfr(:,1:2));
for i=find(iflip>0)';
        im(:,:,i)=circshift(flipud(im(:,:,i)),[mod(n+1,2) 0]);  % flip before rotate.
end;
aliImgs=rsRotateImage(im,-theta);

