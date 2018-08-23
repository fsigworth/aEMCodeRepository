function fvol=OversamplingReconstr2(FourierStack,angles)
% function vol=OversamplingReconstr2(FourierStack,angles)
% Given a stack of nim Fourier transforms of projections (origin in the
% center, n x n pixels each)
% and a 3 x nim matrix of corresponding Euler angles, insert these
% projections into an oversampled 3D Fourier volume.  The volume is then downsampled
% to n^3.  The image size n must be even.



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
x=x(pts);  % select only the active points.
y=y(pts);
ct2=nx/2;  % center offset in oversampled volume

for i=1:nim  % accumulate all the planes
    Im=FourierStack(:,:,i);
    r=EulerMatrix(angles(i,:));  % pick up phi and theta
    X=r(1,1)*x+r(1,2)*y;
    Y=r(2,1)*x+r(2,2)*y;
    Z=r(3,1)*x+r(3,2)*y;
    
    % Compute coords of vol as 1d array, using nearest-neighbor interp.
    P=1+round(ovs*X+ct2)+nx*(round(ovs*Y+ct2)+nx*round(ovs*Z+ct2));
    %     if any(max(P(:)>numel(volx))) || any(min(P(:)<1))
    %         error(['P out of bounds' num2str([max(P(:)) min(P(:))])]);
    %     end;
    
    % Do the accumulation here
    % The fact that the
    % oversampling is >=2 means that, at the insertion of each plane, no point in
    % the volume is affected by more than one plane pixel.  This makes the
    % accumulation code in Matlab very simple.
    volx(P)=volx(P)+Im(pts);
    
    % (Alternate code, also with x and y not selected by pts)
    %     for j=pts
    %         volx(P(j))=volx(P(j))+Im(j);
    %     end;
end;
rvol=fftshift(ifftn(volx)); % go to real space
rdvol=Crop(rvol,n);         % extract the center
fvol=fftn(ifftshift(rdvol)); % back to small Fourier space
