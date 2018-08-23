function out=DownsampleNearest(in,nout,mag,stack)
% function out=DownsampleNearest(in,nout,mag,stack)
% Nearest-neighbor vector, image or volume resampling.
% nout and mag can be scalars, or vectors of size ndim(in).
% If stack is true, 3- or 4-dimensional input is interpreted as a stack of
% images or volumes.  Stacks of 1d vectors is not supported.

% Changed 25 May 15 to maintain centers of images (FFT convention)

nin=size(in);
if isa(in,'numeric')
    inClass=class(in);
else
    inClass='single';  % logical etc. gets convered to single.
end;

% Determine the number of dimensions and extend nout
ndim=numel(nin);

if nargin<4
    stack=false;
end;

if stack  % last dimension will be the number of images
    nImgs=nin(ndim);
    ndim=ndim-1;
    nin=nin(1:ndim);
else
    nImgs=1;
end;

% Convert nout into a vector of the right size
if numel(nout)<ndim
    nout(end:ndim)=max(nout);
elseif numel(nout)>ndim
    nout(ndim+1:end)=[];
end;

% Assign the default mag
if nargin<3 || numel(mag)<1 || any(mag==0)
    mag=nout./nin;
end;

% Extend last elements of mag if needed
if numel(mag)<numel(nin)
    mag(end+1:numel(nin))=mag(end);
end;

if ndim<3  % handle 1d arrays
   nout(nin==1)=1;
end;

% Get the input x corresponding to each output x value
xins=floor(nin(1)/2+1)+round((-floor(nout(1)/2):floor((nout(1)-1)/2))/mag(1));
xins=max(1,min(xins,nin(1)));

yins=floor(nin(2)/2+1)+round((-floor(nout(2)/2):floor((nout(2)-1)/2))/mag(2));
yins=max(1,min(yins,nin(2)));

if ndim<3
    [xq,yq]=ndgrid(xins,yins);
    pins=xq+nin(1)*(yq-1);
    out=zeros([prod(nout) nImgs],inClass);
    in=reshape(in,prod(nin),nImgs);
    for i=1:nImgs
        out(:,i)=in(pins(:),i);
    end;        
    out=reshape(out,[nout nImgs]);
elseif ndim==3
%     3D: assign z's also.
    zins=floor(nin(3)/2+1)+round((-floor(nout(3)/2):floor((nout(3)-1)/2))/mag(3));
    zins=max(1,min(zins,nin(3)));
    [xq,yq,zq]=ndgrid(xins,yins,zins);
    pins=xq+nin(1)*((yq-1)+nin(2)*(zq-1));
    
    out=zeros([prod(nout) nImgs],inClass);
    in=reshape(in,prod(nin),nImgs);
    for i=1:nImgs
        out(:,i)=in(pins(:),i);
    end;
    out=reshape(out,[nout nImgs]);
else
    error(['Dimension of input too large :' num2str(ndim)]);
end;
if isa(in,'logical')
    out=logical(out);
end;
