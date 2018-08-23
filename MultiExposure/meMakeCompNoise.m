function [y cc]=meMakeCompNoise(m1,m2)
% function [y cc]=meMakeCompNoise(m1,m2)
% Given two pre-whitened input images, derive colored noise 
% y to compensate for the fixed pattern
% noise. The raw ccf of m1 and m2 is also returned.
% m1 and m2 should be images that have a large translational shift or
% large defocus difference to avoid interference with true cross
% correlation components.
% y will have an ACF equal to the stored "FixedPatternCCF.mat" but scaled
% by least-squares to the zero-lag CCF of m1 and m2.
% To make m1 and m2 lacking the fixed pattern CCF, add y to m1 and subtract
% y from m2 (or vice versa).  To correct for different doses in new images
% m3 and m4, use y*(d3/d1)+m3 and -y*(d4/d2)+m4

% Pick up the FixedPatternCCF from the CCD library directory.
pa=fileparts(which('CCDMakeFixedPatternRefCCF'));
load([pa '/FixedPatternCCF.mat']);  %loads ccf0
nc=size(ccf0,1);

% Get a cross-correlation function of two images
n=size(m1,1);
win=single(SquareWindow(n,32));
winvar=win(:)'*win(:)/(n^2);

cc=fftshift(real(ifftn(fftn(m1).*win...
                .*conj(fftn(m2).*win))));
ccs=Crop(cc,nc)/winvar;  % normalize by missing window

% Linear least-squares estimation of scaling of the reference ccf
% to match the ccs just computed.
% Set up the matrix eqn ax=y, with x being the coefficient vector
f(:,1)=ccf0(:);
dcf=20*fuzzymask(nc,2,7.1,.01);
f(:,2)=10*dcf(:);  % scale up to make the matrix better conditioned
a=zeros(2,2);
for i=1:2
    for j=1:2
        a(i,j)=f(:,i)'*f(:,j);
    end;
end;
z=ccs(:);
y=f'*z;
coeffs=a\y;
ccf=ccf0*coeffs(1);

% coeffs

% Make some random noise
rn=single(randn(n,n))/n;
fn=fftn(rn);

% Compute the square root of the desired power spectrum.
% Real forces the final CCF to be centrosymmetric.
fc=real(sqrt(fftn(fftshift(Crop(ccf,n)))));

% Create the output function
y=real(ifftn(fn.*fc));
