function S = VesicleAmpGeneral(n,a,s,org)
%  function S = VesicleAmpGeneral(n,a,s,org)
% Like VesicleModelGeneral, but compute the angularly-dependent local image
% signal S as a function in 2D.  We use this in normalizing the particle
% amplitude for particle picking.  We make use of the angular Fourier series
% s (=mi.vesicle.s(ind,:,1)).  The radius a is taken to be a scalar (we
% compute a circular field).  The result is
% S=real(s(1)+s(2)*exp(iTheta)*d1(r)+s(3)*exp(2iTheta)*d2(r)+...)
% where d1(r), d2(r) etc. are cosine-based windows in r.

if numel(n)<2  % scalar means a square image
    n=[n n];
end;
ctr=ceil((n+1)/2);  % correct also for odd n
if nargin<4
    org=ctr;
end;
a=a(1); % only use the first term of radius expansion
s=s(:); % a column vector
nterms=numel(s);

[r2D,theta2D]=Radius(n,org);
n2=numel(r2D);

% Work with 1d variables
r=reshape(r2D,n2,1);  % 1d radius
theta=reshape(theta2D,n2,1);
d0=sqrt(max(0,.5-.5*cos(pi*r/a)));  % cosine window for angular dependence
d0(r>a)=1;
w0=exp(1i*theta);  % complex exponential
wi=ones([n2 nterms],'single');
% Get powers of the complex exponential
for j=2:nterms
    wi(:,j)=wi(:,j-1).*w0(:).*d0;
end;

S=reshape(real(wi*s),n);  % make a 2d array.
