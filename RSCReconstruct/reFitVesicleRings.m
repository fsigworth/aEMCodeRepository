function [fits,rings]=reFitVesicleRings(si,imgs,mask,radii,widths,inds)
% Make a set of orthogonal functions consisting of rings with Gaussian
% radial components.
%

n=size(imgs,1);
nr=numel(radii);
ctr=ceil((n+1)/2);
if nargin<6
    inds=1:size(imgs,3);
end;
nImgs=numel(inds);
fits=zeros(n,n,nImgs,'single');

for j=1:nImgs
    ind=inds(j);
    img=imgs(:,:,ind);
    rves=si.rVesicle(ind);
    yctr=ctr-si.yClick(ind);  % yClick is always positive
    if (yctr+rves<n || yctr-rves>1) && yctr+rves>1  % vesicle is in frame
        rs=radii+rves;
        R=Radius(n,[ctr yctr]);
        
        rings=zeros(n,n,nr,'single');
        ringsM=zeros(n,n,nr,'single');
        ik=1;
        for i=1:nr
            ring=exp(-((R-rs(i)).^2/(2*widths(i)^2)));
            rm=ring.*mask;
            if sum(rm(:))>.1*n*widths(i)*sqrt(2*pi) % at least 10% is visible
                rings(:,:,ik)=ring;
                ringsM(:,:,ik)=rm;
                ik=ik+1;
            else
% disp([j i])
% whos ringsM
                rings(:,:,ik)=[];
                ringsM(:,:,ik)=[];  % don't fit this ring
            end;
        end;
        nr1=size(ringsM,3);
        if nr1>0
            F=reshape(rings,n^2,nr1);
            FM=reshape(ringsM,n^2,nr1);
            b=LinLeastSquares(FM,img.*mask);
            if ~any(isnan(b)) && nr1==nr
            end;
            fits(:,:,j)=reshape(F*b,n,n);
        end;
    end;
end;

% r0=reshape(rings,n^2,nr);
% r1=qr(r0);
% r2=reshape(r1,n,n,nr);
