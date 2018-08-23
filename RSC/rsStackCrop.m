function [si,imgsOut]=rsStackCrop(si, imgs, nOut)
% Takes a StackInfo structure and stack of images. Crop the images and 
% correspondingly downsample the ctfs.
isStackImgs=~ismatrix(imgs);
isStackCtfs=~ismatrix(si.ctfs);

imgsOut=Crop(imgs,nOut,isStackImgs);
si.ctfs=DownsampleGeneral(si.ctfs,nOut,[],isStackCtfs);

