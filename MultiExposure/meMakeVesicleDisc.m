function m=meMakeVesicleDisc(n,a,width,org);
%  Make a fuzzymask-like disc having the shape of a vesicle. The edge of
%  the disc corresponds to the membrane center. The width is the 10-90%
%  risetime of the edge in pixels. If width=0, we return a binary image.
%  The disc approaches the value 1 in the center, 0 outside.


if nargin<3
    width=real(a(1))/10;
end;
if nargin<4
    org=floor((n+1)/2)*[1 1];
end;
if numel(n)<2
    n=n*[1 1];
end;
nTerms=find(abs(a)>0,1,'last');
a=a(1:nTerms);
a=a(:);

k=1.782./width;

[r2D,theta2D]=Radius(n,org);

% Work with 1d variables
n2=prod(n);
r=reshape(r2D,n2,1);  % 1d radius
theta=reshape(theta2D,n2,1);
w=exp(1i*theta);  % complex exponential
wi=ones([n2 nTerms],'single');
% Get powers of the complex exponential
for j=2:nTerms
    wi(:,j)=wi(:,j-1).*w(:);
end;

rNom=real(wi*a);  % nominal radius the shell, 2D function, as a function of theta

if width>0
    m=0.5*(1-erf(k*(r-rNom)));
else
    m=(r-rNom)>=0;
end;
m=reshape(m,n);
