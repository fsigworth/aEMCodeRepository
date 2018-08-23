function gP=gFourierShift(n,sh)
% function gP=gFourierShift(n,sh)
% Compute the complex exponential for shifting an image.  n is the size of
% the image (scalar for a square image or a 2-element vector for a
% rectangular image), and can be odd.  If n is a 3-element vector, then a
% 3D array is produced for shifting a volume. sh =[dx dy] or [dx dy dz] contains the shifts in
% pixels. Positive values are shifts up and to the right.  In the returned
% matrix, P(1,1) is zero frequency.
% If sh is nim x 2 in size and numel(n)<3 then a stack of 2D complex images is created.
%
% e.g. to shift by dx, dy pixels,
%   fm=fftn(m);
%   P=FourierShift(size(m),[dx dy]);
%   fsh=real(ifftn(fm.*P));
% or, in one line,
%   fsh=real(ifftn(fftn(m).*FourierShift(size(m,1),[dx dy])));
%
n=single(n);
sh=single(sh);
if numel(n)<2
    n(2)=n(1);
end;
if numel(sh)<4  % Too small to be a stack of shifts
    sh=sh(:)';  % force a row vector;
end;

ndims=numel(n);
nim=size(sh,1); % one row for each output image
switch ndims
    case 2
        gP=zeros([n nim],'single','gpuArray');
        for i=1:nim
            nx=gpuArray(n(1));
            ny=gpuArray(n(2));
            [gX,gY]=ndgrid((-floor(nx/2):floor((nx-1)/2))/nx,...
                (-floor(ny/2):floor((ny-1)/2))/ny);
            gP(:,:,i)=exp((-1j*2*pi)*fftshift((sh(i,1)*gX+sh(i,2)*gY)));
        end;
    case 3
        if nim>1
            error('Can''t handle a stack of 3D volumes');
        end;
        nx=gpuArray(n(1));
        ny=gpuArray(n(2));
        nz=gpuArray(n(3));
        [gX,gY,gZ]=ndgrid((-floor(nx/2):floor((nx-1)/2))/nx,...
            (-floor(ny/2):floor((ny-1)/2))/ny,...
            (-floor(nz/2):floor((nz-1)/2))/nz);
        gP=exp((-1j*2*pi)*fftshift((sh(1)*gX+sh(2)*gY+sh(3)*gZ)));
    otherwise
        error('n must be a 1, 2 or 3-element vector');
end;

