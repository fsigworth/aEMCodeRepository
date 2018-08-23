function rm = NRotate3(m, rot, origin)
% NROTATE3  Fast map rotation according to the matrix rot,
% using nearest neighbor interpolation.
% rm = Nrotate3(m, rot, origin); rotation about
% origin = [x y z] (default: n/2+1 for x, y and z.  We assume a coordinate system
% with m = m(x,y,z) and x is the fastest-changing index.

[nx, ny, nz]=size(m);
if nargin < 3
    origin=ceil((size(m)+1)/2);  % fft center
end;

N=nx*ny*nz;

% Create input and output coordinate matrices
% First make vectors corresponding to each axis.
rampx= 1-origin(1):nx-origin(1);
rampy= 1-origin(2):ny-origin(2);
rampz= 1-origin(3):nz-origin(3);

% Make 3D cartesian arrays
[x,y,z]=ndgrid(rampx, rampy, rampz);

% Perform the Euler rotation
X=rot(1,1)*x+rot(1,2)*y+rot(1,3)*z+origin(1);
Y=rot(2,1)*x+rot(2,2)*y+rot(2,3)*z+origin(2);
Z=rot(3,1)*x+rot(3,2)*y+rot(3,3)*z+origin(3);
% x=[]; y=[]; z=[];

% Convert to 1D arrays
X=reshape(X,N,1);
Y=reshape(Y,N,1);
Z=reshape(Z,N,1);

% Clip the coordinates to be in bounds.
% eps=1e-6; % We insure that rounding up and down leaves us in bounds.

X=max(X,1);
X=min(X,nx);

Y=max(Y,1);
Y=min(Y,ny);

Z=max(Z,1);
Z=min(Z,nz);

X0=round(X);
% Xi=X-X0;
% X1=X0+1;

Y0=round(Y);
% Yi=Y-Y0;
% Y1=Y0+1;
Y0=(Y0-1)*nx;  % Convert to linear index component for rows.
% Y1=(Y1-1)*nx;

Z0=round(Z);
% Zi=Z-Z0;
% Z1=Z0+1;
Z0=(Z0-1)*nx*ny;  % Convert to linear index component for planes.
% Z1=(Z1-1)*nx*ny;


% No interpolation:
rm= m(X0+Y0+Z0);

% % Linear interpolation:
% rm= (1-Zi).*((1-Yi).*((1-Xi).*m(X0+Y0+Z0) + Xi.*m(X1+Y0+Z0))...
% 				+Yi.*((1-Xi).*m(X0+Y1+Z0) + Xi.*m(X1+Y1+Z0)))...
% 	  + Zi.*((1-Yi).*((1-Xi).*m(X0+Y0+Z1) + Xi.*m(X1+Y0+Z1))...
% 				+Yi.*((1-Xi).*m(X0+Y1+Z1) + Xi.*m(X1+Y1+Z1)));

rm=reshape(rm,nx,ny,nz);
