function [val, coords, A]=max3di(m,shiftToOrigin)
% function [val, coords, A]=max3di(m)
% function [val, coords, A]=max3di(m, shiftToOrigin)
% find the interpolated maximum value and its interpolated
% maximum location in the matrix m(i,j,k).  Interpolation is quadratic and
% uses InterpMax3.
% If shiftToOrigin (default 0) is set to 1, the returned 3x1 vector coords
% will be zero if the peak of m is in the center.

    shiftToOrigin=0;

if nargin<2
else
    shiftToOrigin=floor((size(m)/2+1));
end;

[nx ny nz]=size(m);
[v0 i0 j0 k0]=max3d(m);

i0=max(i0-1,1);
i0=min(i0,nx-2);
j0=max(j0-1,1);
j0=min(j0,ny-2);
k0=max(k0-1,1);
k0=min(k0,nz-2);

[val ctr A]=InterpMax3(m(i0:i0+2,j0:j0+2,k0:k0+2));
coords=[i0-1 j0-1 k0-1]+ctr-shiftToOrigin;
