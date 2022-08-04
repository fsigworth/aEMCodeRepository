function r=radius(n,center)
% returns an nxn matrix containg the distance from the center.
% By default, center=[(n+1)/2 (n+1)/2]
%
if nargin<2
    center = [(n+1)/2 (n+1)/2];  % Coordinate of the center
end;

[x,y]=ndgrid(-center(1)+1:n-center(1),-center(2)+1:n-center(2));
r=sqrt(x.^2+y.^2);
