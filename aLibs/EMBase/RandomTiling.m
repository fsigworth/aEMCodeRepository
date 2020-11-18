function [xs,ys,nCollisions]=RandomTiling(nm,mSize,minDist,mBorder);
% Return nm random positions (ix,iy) inside a rectangle of mSize,
% not allowing positions to be closer than minDist, and
% (optionally) not closer than mBorder to the edges.
maxOver=100; % maximum oversampling allowed.
if nargin<4
    mBorder=0;
end;
areaFraction=nm*minDist.^2/prod(mSize-2*mBorder);
if areaFraction>.5
    warning(['Area fraction is high: ' num2str(areaFraction)]);
end;
ind=1;
nCollisions=0;
xs=zeros(nm,1,'single');
ys=zeros(nm,1,'single');
while ind<nm
    tx=mBorder+(mSize(1)-2*mBorder)*rand;
    ty=mBorder+(mSize(2)-2*mBorder)*rand;
    if ind>0
        dx=xs(1:ind)-tx;
        dy=ys(1:ind)-ty;
        dMin=sqrt(min(dx.^2+dy.^2));
    else
        dMin=inf;
    end;
    if dMin>minDist
        xs(ind,1)=tx;
        ys(ind,1)=ty;
        ind=ind+1;
    else
        nCollisions=nCollisions+1;
        if nCollisions>nm*maxOver
            error(['Too many collisions: ' num2str(nCollisions)]);
    end;
    end;
end;
if areaFraction>.5 % guve a value of nCollisions
        nCollisions
end;
