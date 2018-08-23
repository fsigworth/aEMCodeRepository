function eigenImages=spMakeEigenimages(stack,nFactors,stride,mask)
% function eigenImages=spMakeEigenimages(stack,nFactors,stride,mask)
% From the stack of images, construct nFactors eigenimages.
% Typical use: eigs=spMakeEigenimages(stack,nFactors);
% The stride and mask variables are optional. For very large datasets you
% can set stride to a number >1 to skip images in the stack.  You can also
% give a mask for the images (same size as an image); the default is a
% fuzzy disc with radius 0.4 x image width.

[nx, ny, nim]=size(stack);
if nargin<3
    stride=1;
end;
if nargin<4
    mask=fuzzymask([nx ny],2,.4*nx,.1*nx);
end;
nimgs=1+floor((nim-1)/stride);
maskThreshold=0.1;

pix=mask>maskThreshold;  % We will extract just these pixels.

imgs1d=zeros(sum(pix(:)),nimgs);  % Images as 1D vectors
for i=1:nimgs
    im=stack(:,:,1+(i-1)*stride).*mask;
    imgs1d(:,i)=double(im(pix));  % each image is a column
end;

tic
if nFactors <38
    disp('svds...');
    [u, s, v]=svds(imgs1d,nFactors);
else
    disp('svd...');
    [u, s, v]=svd(imgs1d);
end;
disp('done.');
toc
% u=u(:,1:nfactors);
% s=s(1:nfactors,1:nfactors);
% v=v(:,1:nfactors);
% u has eigenimages as columns, nx*ny by ndims
% s*v' is nfactors by nimgs, the vectors.
% factors=s*v';
eigenImages=zeros(nx,ny,nFactors,'single');
eig1=zeros(nx,ny);
for i=1:nFactors  % Restore the pixels to images.
    eig1(pix)=u(:,i);
    eigenImages(:,:,i)=eig1;
end;
