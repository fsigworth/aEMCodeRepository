function W = VesicleMaskGeneral(n,a,z,org)
% function W = VesicleMaskGeneral(n,a,z,org)
% Given a cross-section coordinate z (0 is center of the membrane), 
% construct a binary mask (1 inside, 0 outside) of a
%  vesicle in an image of size n, with radius a and center org.  All
%  coordinates are in pixels.
% Unlike VesicleFromModel, this function accepts a as a vector, so that
% a(1) is the radius, re(a(2)) is the perturbation in radius as cos(theta),
% re(a(2)) is the sin(theta) perturbation; a(3) is the cos(2theta)
% perturbation (ellipse) and so on.


% Test code
% n=256;
% q=1i;
% a=[80 0 0 10*q 1*q 20 2];
% a=[40 0 0 q q 20 2];
% % a=[40 50 0 0 0 0 0];
% z=0;
% nargin=3;

if numel(n)<2  % scalar means a square image
    n=[n n];
end;
ctr=ceil((n+1)/2);  % correct also for odd n
if nargin<4
    org=ctr;
end;
a=a(:); % make it a column vector;
nterms=numel(a);

% Handle the case that it is a spherical vesicle after all
% --no high-order terms in radius
if nterms<2 || all(a(2:end)==0)
    W=SphereDensity(n,a(1)+z,org)>0;
    return
end;

% Compute n1, the size of the working region
% aExt=1.2*a;  % increased values for computing n1
% aExt(1)=a(1);
org=org(:)';  % make it a row vector
% n1=[1 1]*2*ceil(max((sum(abs(aExt))))+numel(model)*.7)+2;  % minimim square that bounds the vesicle.
% % % if n1<min(n)  % We'll model the vesicle in a smaller square area
% % %     ctr1=ceil((n1+1)/2);
% % %     fracShift=org-round(org);
% % %     org1=fracShift+ctr1;
% % % else
% % %     org1=org;
% % %     n1=n;
% % % end;
n1=n;
org1=org;

n2=prod(n1);  % linear dimension

% Initialization code from SphereDensity
[r2D,theta2D]=Radius(n1,org1);

% Work with 1d variables
r=reshape(r2D,n2,1);  % 1d radius
theta=reshape(theta2D,n2,1);
w=exp(1i*theta);  % complex exponential
wi=ones([n2 nterms],'single');
% Get powers of the complex exponential
for j=2:nterms
    wi(:,j)=wi(:,j-1).*w(:);
end;

rNom=real(wi*a);  % nominal radius of shells, 2D function, as a function of theta

% Correction of membrne thickness for radius variations
rDeriv=real(wi*(1i*a.*(0:nterms-1)'));  % derivative of a wrt theta
thkFactor=sqrt((rDeriv./rNom).^2+1);


    r0=rNom+z*thkFactor;  % corrected radius as a fcn of theta.
        W1=r<r0;
W=reshape(W1,n1);  % convert back to 2D

% imags(W);