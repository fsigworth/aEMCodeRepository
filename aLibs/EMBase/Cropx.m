function mc=Cropx(m,n,isstack,fillval)
% function mc=Cropx(m,n,isstack,fillvalue)
% Reduce the size of the 1d array, square or cube m by cropping
% (or increase the size by padding with fillval, by default zero)
% to a final size of n x n or n x n x n.  Unlike Crop() we do no shift, so
% m(1,1) remains unchanged; instead we add or remove values at the ends.
% If m is 2-dimensional and n is a vector, m is cropped to n=[nx ny].
% If the flag isstack = 1 then a 3D array m is treated as a stack of 2D
% images, and each image is cropped to n x n.
% For 2D images, the input image doesn't have to be square.
% The result is double if fillval is double; by default the result is
% single.
% Now handles simultaneous cropping/padding for 2D input. fs Apr-20.
% 3D-4D inputs are still
% assumed to be stacks of square, or cubic in dimension.

if nargin<3
    isstack=0;
end;
if nargin<4
    fillval=single(0);  % Force a single output when padding.
else
    fillval=single(fillval);
end;

szm=size(m);
ndi=ndims(m);
if ndi==2 && any(szm==1) % 1D array
    ndi=1;
    n=[max(n) 1];
    m=reshape(m,n); % force a column vector
else
    if numel(n)<ndi % extend the explicit number of dimensions
        n(end+1:ndi)= n(end);
    end;
end;

if isstack
    nim=szm(ndi); % copy all the slices
    ndi=ndi-1;
    n=n(1:ndi);
    szm=szm(1:ndi);
else
    nim=1;
end;

no=min(n,szm); % number of points to copy

if any(n>szm) % any padding
    mc=zeros([n nim],'single');
end;
for iSlice=1:nim
    switch ndi
        case 1
            mc(1:no,iSlice)=m(1:no,iSlice);
        case 2
            mc(1:no(1),1:no(2),iSlice)=m(1:no(1),1:no(2),iSlice);

        case 3 % m is 3D
            mc(1:no(1),1:no(2),1:no(3),iSlice)=m(1:no(1),1:no(2),1:no(3),iSlice);
        otherwise
            error(['Cropx: dimension too large: ' num2str(ndi)]);
    end;
end;
