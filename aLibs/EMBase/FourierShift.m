function P=FourierShift(n,sh)
% function P=FourierShift(n,sh) Compute the complex exponential for
% shifting an image.  n is the size of the column vector, image or volume
% (can be a scalar for a square image or a cube), and can be odd. sh is the
% desired shift in pixels, and its dimension sets the dimension of P.
% Positive values are shifts up and to the right.  In the returned
% matrix, P(1) is zero frequency. If sh is nim x ndims in size then a
% stack of vectors/images/volumes is created.
%
% e.g. to shift an image by by dx, dy pixels,
%   fm=fftn(m);
%   P=FourierShift(size(m),[dx dy]);
%   fsh=real(ifftn(fm.*P));
% or, in one line,
%   fsh=real(ifftn(fftn(m).*FourierShift(size(m,1),[dx dy])));
%

nim=size(sh,1); % multiple rows for a stack of outputs
ndims=size(sh,2); % columns of shifts set the dimension
if numel(n)<ndims
    n=ones(1,ndims)*n(1);
end;
        P=zeros([n nim]);
switch ndims
    case 1
       X=((-floor(n/2):floor((n-1)/2))/n)';
for i=1:nim
            P(:,:,i)=exp((-1j*2*pi)*fftshift(sh(i,1)*X));
end;
    case 2
        for i=1:nim
            nx=n(1);
            ny=n(2);
            [X,Y]=ndgrid((-floor(nx/2):floor((nx-1)/2))/nx,...
                (-floor(ny/2):floor((ny-1)/2))/ny);
            P(:,:,i)=exp((-1j*2*pi)*fftshift((sh(i,1)*X+sh(i,2)*Y)));
        end;
    case 3
        if nim>1
            error('Can''t handle a stack of 3D volumes');
        end;
        nx=n(1);
        ny=n(2);
        nz=n(3);
        [X,Y,Z]=ndgrid((-floor(nx/2):floor((nx-1)/2))/nx,...
            (-floor(ny/2):floor((ny-1)/2))/ny,...
            (-floor(nz/2):floor((nz-1)/2))/nz);
        P=exp((-1j*2*pi)*fftshift((sh(1)*X+sh(2)*Y+sh(3)*Z)));
    otherwise
        error('n must be a 1, 2 or 3-element vector');
end;
