function h=imaga(X,Y,m)
% function imaga(m);
% Draw an image using cartesian coordinates; the image pixel values are 
% assumed to lie in 0...255.It is rendered in 256-level
% grayscale regardless of the colormap in use.  The coordinates are
% cartesian with m(1,1) being in the lower left, m(end,1) on the lower right.
% If the arguments X and Y are given, these define the
% x and y axes as in image().

if nargin<2    
    z=squeeze(X)';
    h=image(repmat(single(z)/256,1,1,3)); axis xy;
else
    z=squeeze(m)';
    h=image(X,Y,repmat(single(z)/256,1,1,3)); axis xy;
end;
if nargout<1
    clear h  % Don't echo the handle if no output variable given.
end;