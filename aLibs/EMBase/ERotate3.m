function rm = ERotate3(m, angles, origin, nearest)
% EROTATE3  Fast map rotation through Euler angles.
% rm = Erotate3(m, angles, origin, nearest): rotate the 3d map m according to the
% rm = ERotate3(m,A,[],nearest) % apply the affine transform A. A is 4x4 or (with no
% In the first form, given are Euler angles and the origin about which
% rotations are done. Alternatively, A is taken to be a 3x3 rotation matrix
% or a 4x4 affine matrix), in which case the origin is taken from the last
% column.
%
% Euler angles are defined as angles = [theta phi psi], about
% origin = [x y z] (default: ceil(n/2+1), the FFT center, for x, y and z.
% We assume a coordinate system
% with m = m(x,y,z) and x is the fastest-changing index.
% By default (nearest = 0) trilinear interpolation is used. If nearest=1,
% nearest-neighbor interpolation is used.

[nx,ny,nz]=size(m);

if nargin<3 || numel(origin)==0
    origin=ceil((size(m)+1)/2);  % fft center
end;
if nargin<4
    nearest=0; % linear interpolation is the default
end;

% The second argument can be a vector of angles, a 3x3 rotation matrix, or
% a 4x4 affine matrix.
outShift=[0 0 0];
sza=size(angles);
switch prod(sza)
    case 3 % a vector
        rot=EulerMatrix(angles(:));
    case 9 % presumably 3x3
        if ~all(sza==[3 3])
           error(['Rotation matrix has the wrong shape: ' num2str(sza)]);
        end;
        rot=angles;
    case 16
        if ~all(sza==[4 4])
           error(['Affine transform matrix has the wrong shape: ' num2str(sza)]);
        end;
        rot=angles(1:3,1:3);
%         origin=origin(:);
        outShift=-rot\angles(1:3,4); % last column is the
%         relative shift added after rotation.
    otherwise
           error(['Angles / transform matrix has the wrong size: ' num2str(sza)]);
           
end;
rot=inv(rot);
% Create input and output coordinate matrices
% First make vectors corresponding to each axis.
rampx= 1-origin(1):nx-origin(1); % zero at index = origin
rampy= 1-origin(2):ny-origin(2);
rampz= 1-origin(3):nz-origin(3);

% Make 3D cartesian arrays
[x,y,z]=ndgrid(rampx, rampy, rampz);

% Perform the Euler rotation
X=rot(1,1)*x+rot(1,2)*y+rot(1,3)*z+origin(1)+outShift(1);
Y=rot(2,1)*x+rot(2,2)*y+rot(2,3)*z+origin(2)+outShift(2);
Z=rot(3,1)*x+rot(3,2)*y+rot(3,3)*z+origin(3)+outShift(3);

% Convert to 1D and clip the coordinates to be in bounds.
epsl=1e-6; % We insure that rounding up and down leaves us in bounds.
epsh=1+1e-6;
X=max(min(X(:),nx-epsh),1+epsl);
Y=max(min(Y(:),ny-epsh),1+epsl);
Z=max(min(Z(:),nz-epsh),1+epsl);

if nearest
    % No interpolation:
    rm= m(round(X)+nx*(round(Y-1)+ny*round(Z-1)));
else
    % Linear interpolation:
    X0=floor(X);
    Xi=X-X0;
    X1=X0+1;
    
    Y0=floor(Y);
    Yi=Y-Y0;
    Y1=Y0+1;
    Y0=(Y0-1)*nx;  % Convert to linear index component for rows.
    Y1=(Y1-1)*nx;
    
    Z0=floor(Z);
    Zi=Z-Z0;
    Z1=Z0+1;
    Z0=(Z0-1)*nx*ny;  % Convert to linear index component for planes.
    Z1=(Z1-1)*nx*ny;
    
    rm= (1-Zi).*((1-Yi).*((1-Xi).*m(X0+Y0+Z0) + Xi.*m(X1+Y0+Z0))...
        + Yi.*((1-Xi).*m(X0+Y1+Z0) + Xi.*m(X1+Y1+Z0)))...
        +Zi.*((1-Yi).*((1-Xi).*m(X0+Y0+Z1) + Xi.*m(X1+Y0+Z1))...
        + Yi.*((1-Xi).*m(X0+Y1+Z1) + Xi.*m(X1+Y1+Z1)));
end;
rm=reshape(rm,nx,ny,nz);
