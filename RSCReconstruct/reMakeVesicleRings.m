function [residuals,a,rings]=reFitVesicleRings(si,imgs,mask,radii,widths,inds)
% Make a set of orthogonal functions consisting of rings with Gaussian
% radial components.
%
n=size(imgs,1);
nr=numel(radii);
ctr=ceil((n+1)/2);
if nargin<4
    inds=1:size(imgs,3);
end;
nImgs=numel(inds);
residuals=zeros(n,n,nImgs,'single');
a=zeros(nr,nImgs);

for j=1:nImgs
    ind=inds(j);
    rs=radii+si.rVesicle(ind);
    yctr=ctr-si.yClick(ind);
    R=Radius(n,[ctr yctr]);
    
    rings=zeros(n,n,nr,'single');
    for i=1:nr
        rings(:,:,i)=mask.*exp(-((R-rs(i)).^2/(2*widths(i)^2)));
    end;
    
    F=reshape(rings,n^2,nr);
    a(:,j)=LinLeastSquares(F,imgs);
    residuals(:,:,j)=img-reshape(F*a(:,j),n,n);
end;

% r0=reshape(rings,n^2,nr);
% r1=qr(r0);
% r2=reshape(r1,n,n,nr);
