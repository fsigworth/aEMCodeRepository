function slices=OversamplingProject2(vol,angles,ovs)
% function slices=OversamplingProject(vol,angles,ovs)
% Extract an n x n 2D slice oriented according to the Euler angles, from a 3D
% Fourier volume of size 4n x 4n x 4n.  angles is nim x 3 in size.
% Nearest-neighbor 'interpolation' is used.  ovs=4 by default.

if nargin<3
    ovs=4;
end;

% Test code
% n=64;
% nx=ovs*n;
% rvol=fuzzymask(n,3,n*0.05,n*0.02,[n*.4 n/2+1 n/2+1]);
% rxvol=Crop(rvol,nx);
% vol=fftshift(fftn(rxvol));
% angles=[0 0 0];
%
nproj=size(angles,1);
nx=size(vol,1);  % oversampled volume
n=nx/ovs;        % size of the returned slice, which is normally sampled
[x y]=ndgrid(-n/2:n/2-1);  % coordinates of the slice plane
pts=find(x.^2+y.^2<(n/2-1)^2);  % valid points
x=x(pts);  % select only the valid points
y=y(pts);

ctx=nx/2;  % offset for center of higher dimensions

slice=complex(zeros(n,n));
slices=complex(zeros(n,n,nproj));

for i=1:nproj
    r=EulerMatrix(angles(i,:)); % Do the rotations
    X=r(1,1)*x+r(1,2)*y;
    Y=r(2,1)*x+r(2,2)*y;
    Z=r(3,1)*x+r(3,2)*y;
    
    % Compute coords of vol as 1d array
    P=1+round(ovs*X+ctx)+nx*(round(ovs*Y+ctx)+nx*round(ovs*Z+ctx));
    %     if any(max(P(:)>numel(vol))) || any(min(P(:)<1))
    %         error(['P out of bounds' num2str(max(P(:))) '  ' num2str(min(P(:)))...
    %             '  ' num2str(numel(vol))]);
    %     end;
    %     s=vol(P).*mask;
    slice(pts)=vol(P);
    slices(:,:,i)=slice;
end;
