function [imgs,doses,effCtf,mc]=meMakeMergeImageSet(mi,cpe,ds,extraFilter,nZeros,skipOutliers,mode)
% Load the original images, transform, downsample, and create a set of
% images ready to add together to merge.  The merged image mc is the same as sum(imgs,3).
% extraFilter (default: 1) is the same size as the output images, center zero.
% This could be for example the 
% The corresponding vesicle model is created by meMakeModelVesicleSet.
if nargin<4
    extraFilter=1;
end;
if nargin<5
    nZeros=1;
end;
if nargin<6
    skipOutliers=0;
end;
if nargin<7
    mode=1;
end;
nim=numel(mi.imageFilenames);
fnames=cell(1,nim);
disp('Reading images:');
for j=1:nim
    fnames{j}=[mi.basePath mi.imagePath mi.imageFilenames{j}];
end;

[m,~,doses]=meReadImagesNorm(fnames,cpe,0,mi.pixA,mi.weights,skipOutliers);
modelSpectrum=CCDModelSpectrum2D(mi.camera);  % handle DE-12 or CCD
m=mePreWhiten(m,modelSpectrum);
[nx, ny, nim]=size(m);
n=[nx ny]/ds;
%
[effCtf, mc, mts, coeffs]=meCombineImages2(-m,mi,ds,1,nZeros,mode);
imgs=single(zeros([n nim]));
for i=1:nim
    imgs(:,:,i)=real(ifftn(ifftshift(extraFilter.*coeffs(:,:,i))...
        .*fftn(mts(:,:,i))));
end;
