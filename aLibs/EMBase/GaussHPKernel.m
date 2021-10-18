function H=GaussHPKernel(sz,fc)
% function out=GaussHP(sz,fc)
% GaussHPKernel: Kernel of n-dimensional HP Gaussian filter, of size sz.
% H is returned with zero frequency at center. So to make a HP filter, do
%     out=real(ifftn(fftn(in).*ifftshift(H)));
% fc is the half-power frequency in units of the sampling frequency.  A
% typical value of fc to give 2 cycles across the width w of the image is
% 2/w.

if fc==0  % Zero frequency: don't filter
    H=1;
    return
end;

ndim=sum(sz>1);

k=-(log(2)/2)*fc^2;
m=sz;
if ndim==1
    x=(-1/2:1/n:1/2-1/n)';
    if sz(1)==1 % a row vector
        x=x';
    end;
        H=exp(k./(x.^2+eps));	% Inv Gaussian kernel
elseif ndim==2
    [x,y]=ndgrid(-1/2:1/m(1):1/2-1/m(1), -1/2:1/m(2):1/2-1/m(2));
    H=exp(k./(x.^2+y.^2+eps));	% Inv Gaussian kernel
elseif ndim==3
    [x,y,z]=ndgrid(-1/2:1/m(1):1/2-1/m(1), -1/2:1/m(2):1/2-1/m(2), -1/2:1/m(3):1/2-1/m(3));
    H=exp(k./(x.^2+y.^2+z.^2+eps));	% Inv Gaussian kernel
else 
    error('Size has improper dimension or > 3');
end;
