function h = CTF3(n, res, lambda, defocus, Cs, B, alpha)
% f = CTF3(n, res, lambda, defocus, Cs, B, alpha)
% Compute the 3D contrast transfer function corresponding
% to the given resolution 'res' in �/pixel; Cs is in mm,
% B in �^-2 and alpha in radians.  Defocus is in microns.
% The result is returned in an nxnxn matrix with h(n/2+1) corresponding
% to the zero-frequency amplitude.  Thus you must use fftshift() on
% the result before performing a Fourier transform.  For example,
% to simulate the effect of the CTF on an image m, do this:
% fm=fftshift(fft2(m));
% cm=ifft3(fftshift(fm.*ctf()));
%
% Note: at present, Cs is ignored.
%
% The first zero occurs at lambda*defocus*f0^2=1.
% e.g. when lambda=.025A, defocus=1um, then f0=1/16�.

f0 = 1/(n*res);  % Spatial frequency unit (inverse �)

k=pi*lambda*defocus*1e4*f0^2;

kr=f0^2*B;  % B-factor

n2=n/2;
[mx my mz]=meshgrid(-n2:n2-1);
r2=mx.^2+my.^2+mz.^2;
h=sin(-k*r2-alpha).*exp(-kr*r2);
