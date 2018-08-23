function [img, locs]=MakeRSCImage(n,map,angles,r,ctr,iso)
% Place particles in an image of size n x n according to the given angles
% and vesicle radius.  For right-side out particles, r should be set to
% a - membraneOffset, where a is the vesicle radius.  the returned na x 2
% array locs gives the coordinates of the particles.
if nargin<5
    ctr=[1 1]*ceil((n+1)/2);
end;
na=size(angles,1);
if nargin<6
    iso=0;
end;
if numel(iso)<na
    iso=iso(1)*ones(na,1);
end;
img=zeros(n);
pangles=angles;
for i=1:na
    if iso(i)
        ang=angles(i,:);
        pangles(i,:)=mod([180+ang(1) 180-ang(2) ang(3)],360);
    end;
end;
templates=rsMakeTemplates(pangles,map);
locs=zeros(na,2);
for i=1:na
    alpha=angles(i,1);
    beta=angles(i,2);
    a=r*sind(beta);
    x=a*sind(alpha);
    y=a*cosd(alpha);
    locs(i,:)=round([x y]+ctr);
    temp=ExtractImage(templates(:,:,i),locs(i,:),n,1);
    img=img+temp;
end;
