function rvol=OversamplingReconstr2(FourierStack,angles,symmetry)
% function vol=OversamplingReconstr2(FourierStack,angles,symmetry)
% Given a stack of nim Fourier transforms of projections (origin in the
% center, n x n pixels each)
% and a 3 x nim matrix of corresponding Euler angles, insert these
% projections into an oversampled 3D volume.  The volume is then downsampled
% to n^3, and the real-space backprojected volume is returned.
% The image size n must be even.

% % test code
% n=80;
% nim=300;
% FourierStack=zeros(n,n,nim);
% FourierStack(:,:,1)=fuzzymask(n,2,n/2-2,2);
% angles=zeros(3,nim);
% for i=1:nim
%     FourierStack(:,:,i)=FourierStack(:,:,1);
%     angles(2,i)=0.2*(i-1);
% end;


[n n2 nim]=size(FourierStack);
ovs=4;  % oversampling factor
nx=ovs*n;
volx=zeros(nx,nx,nx);

[x y]=ndgrid(-n/2:n/2-1); % Coordinates in the slice plane
pts=find((x.^2+y.^2)<(n/2-1)^2);  % active points in the plane
x=x(pts);  % select only the active ones.
y=y(pts);
ct2=nx/2;  % center offset in oversampled volume

for i=1:nim  % accumulate all the planes
    img=FourierStack(:,:,i);
    angs=angles(i,:);
    for j=0:symmetry-1
        sangs=angs;
        sangs(1)=angs(1)+2*pi*j/symmetry;

% We take the inverse matrix because we are rotating the plane, not the
% coordinate system.
    r=inv(EulerMatrix(sangs));  % pick up phi and theta
    X=r(1,1)*x+r(1,2)*y;
    Y=r(2,1)*x+r(2,2)*y;
    Z=r(3,1)*x+r(3,2)*y;

    % Compute coords of vol as 1d array, using nearest-neighbor interp.
    P=1+round(ovs*X+ct2)+nx*(round(ovs*Y+ct2)+nx*round(ovs*Z+ct2));
    % Do the accumulation here
    % The fact that the
    % oversampling is 4 means that, at the insertion of each plane, no point in
    % the volume is affected by more than one plane pixel.  This makes the
    % accumulation code in Matlab very simple.
    volx(P)=volx(P)+img(pts);
    end;
end;
rvolx=fftshift(ifftn(ifftshift(volx))); % go to real space
% rvolx=fftshift(ifftn((volx))); % go to real space
rvol=Crop(rvolx,n);         % extract the center
