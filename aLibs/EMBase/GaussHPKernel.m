function ker=GaussHPKernel(n, fc)
% function ker=GaussHPKernel(n, fc)
% Create the kernel (zero frequency in the center) for a Gaussian HP filter.
% The dimension of ker is set by the number of non-singleton elements of n.
% The GaussHP function could be implemented as
%   n=size(in);
%   ker=GaussHPKernel(n,ndims,fc);
%   out=real(ifftn(fftn(in).*ker));
% The kernel is for filtering the input image to give the half-power frequency fc
% (in units of the sampling frequency).  A typical value of fc
% to give 2 cycles across the width n of the image is 2/n.

if fc<.01/max(n) % less than 1% of a cycle across the image width
    ker=ones(n,'single');
else
    k=log(2)/2*fc^2;
    f=RadiusNorm(n);
    ker=exp(-k./f.^2+eps('single'));
end;