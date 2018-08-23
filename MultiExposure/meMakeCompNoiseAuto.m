function [fy cc]=meMakeCompNoiseAuto(m1,m2)
% function [fy cc]=meMakeCompNoiseAuto(m1,m2)
% Given two pre-whitened input images, derive the FT of 
% the colored noise fy to compensate for the fixed pattern
% noise. The raw ccf of m1 and m2 is also returned.
% m1 and m2 should be images that have a large translational shift or
% large defocus difference to avoid interference with true cross
% correlation components.
% Unlike meMakeFixedPatternCompNoise (no Auto) this version does not read
% the stored stndard ACF, but just picks out the center of the computed
% ACF gotten here.
% To make m1 and m2 lacking the fixed pattern CCF, add y to m1 and subtract
% y from m2 (or vice versa).  To correct for different doses in new images
% m3 and m4, use y*(d3/d1)+m3 and -y*(d4/d2)+m4. 
%
% For example, suppose m(:,:,i), i=1..3 consists of three exposures at
% near, middle and far from focus, with doses d(i).
%   fy=meMakeCompNoiseAuto(m(:,:,2),m(:,:,3));
%   d0=sqrt(d(2)*d(3));
%   for i=1:3
%       mcorr(:,:,i)=real(ifftn(fftn(m(:,:,i))+(-1)^i*d(i)/d0*fy));
%   end

nc=32;

% Get a cross-correlation function of two images
n=size(m1);
win=single(SquareWindow(n,32));
winvar=win(:)'*win(:)/prod(n);

cc=fftshift(real(ifftn(fftn(m1.*win)...
                .*conj(fftn(m2.*win)))));
ccs=Crop(cc,nc)/winvar;  % normalize by missing window

dcf=fuzzymask(nc,2,4,1);  % outer mask is radius of 4 pixels
dcc=fuzzymask(nc,2,1.5,.01);  % inner mask is 3x3 square.
annmask=dcf-dcc;
annulus=ccs.*annmask;

% imacs(annmask); drawnow;

ccf=dcf.*(ccs-sum(annulus(:))/sum(annmask(:)));

% Make some random noise
rn=randn(n);
fn=fftn(rn);

% Compute the square root of the desired power spectrum.
% Real forces the final CCF to be centrosymmetric.
fc=real(sqrt(fftn(fftshift(Crop(ccf,n)))/prod(n)));

% Create the output function
% y=real(ifftn(fn.*fc));
fy=fn.*fc;
