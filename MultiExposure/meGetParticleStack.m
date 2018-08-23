function [stack pixA]=meGetParticleStack(m,mi,boxSize,indices)
% function [stack pixA]=meGetParticleStack(m,mi,boxSize,indices)
% Given an image m and its corresponding mi structure, create a stack of
% particle images, extracted according to rounded pixel coordinates.  The
% pixel size of the stack images is computed according to the amount of
% downsampling of m.
% If the optional argument indices is given, only the
% particles with those indices are placed into the stack.  The pixel size of
% the resulting stack is the same as the pixel size of m.  The final stack
% is of size boxSize x boxSize x numel(indices)
% I think it's better not to use the mi.boxSize field, but if only two arguments
% are given, this function uses that field.

if nargin<4
    indices=1:numel(mi.particle.x);
end;
if nargin<3
    boxSize=mi.boxSize;
end

ni=numel(indices);
n=size(m);
ds=mi.imageSize(1)/n(1);  % downsampling factor of the given image
pixA=ds*mi.pixA;

% Stored particle coordinates range from 0..n-1 by convention.
coords=round([mi.particle.x(:) mi.particle.y(:)]/ds+1);
stack=single(zeros(boxSize,boxSize,ni));
for i=1:ni
    p=coords(indices(i),:);
    stack(:,:,i)=ExtractImage(m,p,[boxSize boxSize]);
end;
