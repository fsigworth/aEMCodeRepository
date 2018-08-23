function out=Interpolate1(in,nout,expansion,inctr)
% function out=Interpolate1(in,inctr,nout,expansion)
% 1D interpolation using spline.  Given the input vector in, create an
% output vector with oversampling by the factor (expansion).  the optional argument
% inctr is the index of the input array that is taken to be the center 
% (x=0) value.  the output array has the center at the usual fft position,
% ceil((nout+1)/2).
% For example, suppose we have a membrane profile p sampled at pixA angstroms
% per pixel, and this profile has its center at point 17.  then
% out=Interpolate1(p,80,pixA,17) will give the funciton sampled at 1A per
% sample.

n=numel(in);
if nargin<4
    inctr=ceil((n+1)/2);
end;
octr=ceil((nout+1)/2);
inx=1-inctr:n-inctr;
outx=(1-octr:nout-octr)/expansion;
% Force the output points to be in bounds
outx=min(max(outx,inx(1)),inx(n));
out=interp1(inx,in,outx);
