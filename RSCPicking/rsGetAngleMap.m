function ptrs=rsGetAngleMap(n,r0,iso,inds,ctr)
% Make an n x n array of pointers pointing to the approximate sphere angles
% for each point on the projected hemisphere of radius r0.  Note that r0
% should be either r-membraneOffset (rso) or r+membraneOffset (iso case).
% The angles are modified for inside-out orientation if iso is true.  The
% map is constructed on the basis of the array inds returned by
% rsGetAnglePointers.
% Example:
% % Set up the list of angles, and pointers to these angle values.
%   angles=rsListHemisphereAngles(10,32);
%   % the angles matrix is 2 x 214 in this case.
%   % Look at the angles:
%   x=sin(angles(1,:)).*cos(angles(2,:));
%   y=sin(angles(1,:)).*sin(angles(2,:));
%   plot(x,y,'k.-');
%
%   % Make a look-up for the angles for each point in a projection
%   ptrs=rsGetAngleMap(128,60,0,inds);
%   % For each i,j in the plane, the corresponding angles are
%   % angles(:,ptrs(i,j))
if nargin<6
    ctr=ceil((n+1)/2)*[1 1];
end;

[nAlpha, nBeta]=size(inds);
[x, y]=ndgrid(1-ctr(1):n-ctr(1),1-ctr(2):n-ctr(2));

r=sqrt(x.^2+y.^2);
alpha=90-180/pi*atan2(y,x);
alpha=mod(alpha,360);

b=min(r/r0,1);
beta=real(180/pi*asin(b));
if iso
%     angles=[mod(180+alpha(:),360) 180-beta(:)];
    angles=[mod(180+alpha(:),360) beta(:)];  % we don't modify beta
else
    angles=[alpha(:) beta(:)];
end;
ptrs=rsGetAnglePointers(angles,nAlpha,nBeta,inds);
ptrs=reshape(ptrs,n,n);
