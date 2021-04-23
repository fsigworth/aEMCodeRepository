% function imagsarr(imgs);
% ni=size(imgs,3);
ni=20

imgs=randn(64,64,20);
minv=min(imgs(:));
% imgs=imgs-minv; % shift so zero will be black
n=size(imgs,1);
nr=ceil(sqrt(ni))
nc=floor(ni/nr)
imgs(:,:,nr*nc)=0; % pad
bigImg=minv*ones((n+1)*nc+1,(n+1)*nr+1,'single');
for i=1:ni
    ic=mod(i-1,nc)+1;
    ir=ceil((i-ic)/nr)+1;
    disp([ic ir]);
    shiftx=(n+1)*(ic-1); % start at the top left
    shifty=(n+1)*(nr-ir);
    bigImg(shiftx+2:shiftx+n+1,shifty+2:shifty+n+1)=imgs(:,:,1);
%     ic=i-(nc-1)*ir;
end;
% bigImg=circshift(bigImg,[1,0]);
imags(bigImg);
