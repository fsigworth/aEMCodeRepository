function bigImage=imagsarray(imgs,nr);
% Low-tech version of imagsar for a tiled display of a stack of square
% images. There is no resizing or mouse selection. The optional argument nr
% is the number of rows in the display, in case we aren't computing it
% automatically. If no output argument is given, we go ahead and call
% imags(bigImage).

% imgs=randn(64,64,20);
ni=size(imgs,3);
minv=min(imgs(:));
n=size(imgs,1);
if nargin<2
    nr=ceil(sqrt(ni));
end;
nc=ceil(ni/nr);
% imgs(:,:,nr*nc)=0; % pad
    % When we create the big image, we make the empty space black.
bigImage=minv*ones((n+1)*nc+2,(n+1)*nr+2,'single');
for i=1:ni
    ic=mod(i-1,nc)+1;
    ir=ceil((i-ic)/nc)+1;
    %      disp([ic ir]);
    shiftx=(n+1)*(ic-1); % start at the top left, skip one pixel each time.
    shifty=(n+1)*(nr-ir);
    bigImage(shiftx+2:shiftx+n+1,shifty+2:shifty+n+1)=imgs(:,:,i);
end;
if nargout<1
    imags(bigImage);
end;
