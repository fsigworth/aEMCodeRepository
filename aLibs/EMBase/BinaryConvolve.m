function out=BinaryConvolve(in,mask)
% function out=BinaryConvolve(in,mask,stack)
% binary convolution of 1d or 2d or 2d stack of images
% (expand by or'ing with the mask)
% !!Actually, at present it computes the correlation of in * mask !!!
% form in * mask for logical variables.  
% It's not clear why, but at present mask must be the same size as in.
% It is assumed to be zero-centered with a 
% compact region of 1s near the center.  Out will have regions of 1s enlarged.

% % % test code
% n=64;
% nim1=2;
% nim2=5;
% in=false(n,n,nim1,nim2);
% for i=1:nim1
%     for j=1:nim2;
%     in(:,:,i,j)=GaussFilt(randn(n),.1)>.1;
% end;
% end;
% % mask=fuzzymask(n,2,5,1,[n/2-3 n/2+1])>.5;
% mask=true(3);
% subplot(2,1,1);
% imovie(in);
% 
% % Trying out the built-in function
% % mask=fuzzymask(n,2,5,1,[5 5])>.5;
% % out+conv2(single(in),single((mask)),'same')>0;
% % subplot(2,1,2);
% % imacs(out+in);
% % return

% if nargin<3
%     stack=false;
% end;
sz0=size(in);
sz=sz0(1:2);
if numel(sz0)>2
    nim=prod(sz0(3:end));
else
    nim=1;
end;
in=reshape(in,[sz nim]);
szm=size(mask);
mctr=ceil((szm+1)/2);
mskx=any(mask,2);
minx=find(mskx,1,'first');
maxx=find(mskx,1,'last');

msky=any(mask,1);
miny=find(msky,1,'first');
maxy=find(msky,1,'last');

lowerPads=max(mctr-[minx miny],0);
upperPads=max([maxx maxy]-mctr,0);
padIn=false([sz+lowerPads+upperPads nim]);
padIn(1+lowerPads(1):sz(1)+lowerPads(1),(1+lowerPads(2):sz(2)+lowerPads(2)),:)=in;
nmsk=[maxx-minx maxy-miny]+1;
msk=mask(minx:maxx,miny:maxy);
out=false([sz nim]);
for i=1:nmsk(1)
    for j=1:nmsk(2)
        if msk(i,j)
            out=out|padIn(i:i+sz(1)-1,j:j+sz(2)-1,:);
        end;
    end;
end;
out=reshape(out,sz0);
% % % test code
% subplot(2,1,2);
% in=reshape(in,sz0);
% imovie(out+in,.2);