function W = gpuVesicleFromModelGeneral(n,a,model,org,aw)
% function W = VesicleFromModelGeneral(n,a,model,org,aw)
% Given a density cross-section vector 'model', construct the density of a
%  vesicle in an image of size n, with radius a and center org.
% Unlike VesicleFromModel, this function accepts a as a vector, so that
% a(1) is the radius, re(a(2)) is the perturbation in radius as cos(theta),
% re(a(2)) is the sin(theta) perturbation; a(3) is the cos(2theta)
% perturbation (ellipse) and so on.
% The optional parameter aw is a set of weights multiplying the resulting function,
% as coefficients of a Fourier expansion in theta.  Thus
% W=W0*re(aw(1)+aw(2)*exp(iTheta)+aw(3)*exp(2iTheta)+...) up to the same
% number of terms as in a.

% Test code
% for iq=1:360;
%     q=exp(1i*iq*pi/180);
% n=256;
% a=[80 0 0 10*q 1*q 20 2];
% a=[80 0 0 q q 10 2]/8;
% a=[10 0 0 0 0 0 0];
% model=[1 1 .7 .5 .7 1 1];
% model=[ .7 .5 .7 ];
% aw=[1 0 0 0 0 0 0];


tol = 1e-3;  % minimum density different to invoke calculation of a shell
rBoundFactor=.9; % radius for damping the angular dependence of radius and
%       and amplitude (should be 0 for no damping, 1 for max damping).

if numel(n)<2  % scalar means a square image
    n=[n n];
end;
ctr=ceil((n+1)/2);  % correct also for odd n
if nargin<4
    org=ctr;
end;
if nargin<5
    aw=1;
end;
a=a(:); % make it a column vector;
nterms=numel(a);

if numel(aw)<nterms
    aw(nterms)=0;  % extend the array.
elseif numel(aw)>nterms
    aw=aw(1:nterms);
end;
aw=gpuArray(aw(:));

% Handle the case that it is a spherical vesicle after all
% --no high-order terms in radius or weights
if nterms<2 || (all(a(2:end)==0) && all(aw(2:end)==0))
    W=aw(1)*VesicleFromModel(n,a(1),model,org);
    return
end;

% Compute n1, the size of the working region
aExt=1.2*a;  % increased values for computing n1
aExt(1)=a(1);
org=org(:)';  % make it a row vector
n1=[1 1]*2*ceil(max((sum(abs(aExt))))+numel(model)*.7)+2;  % minimim square that bounds the vesicle.
if n1<min(n)  % We'll model the vesicle in a smaller square area
    ctr1=ceil((n1+1)/2);
    fracShift=org-round(org);
    org1=fracShift+ctr1;
else
    org1=org;
    n1=n;
end;

N1=gpuArray(n1);

N2=prod(N1);  % linear dimension

model=gpuArray([0;model(:);0]);  % pad with zeros
nm=numel(model);
maxDensity=max(abs(model));
mctr=ceil((nm+1)/2);

% Initialization code from SphereDensity
[r2D,theta2D]=Radius(N1,org1);

% Work with 1d variables
r=reshape(r2D,N2,1);  % 1d radius
theta=reshape(theta2D,N2,1);
w=exp(1i*theta);  % complex exponential
wi=ones([N2 nterms],'single','gpuArray');
% Get powers of the complex exponential
for j=2:nterms
    wi(:,j)=wi(:,j-1).*w(:);
end;

% limits of integration over 1 pixel
rm=r-.5;
rm2=rm.^2;
rp=r+.5;
rp2=rp.^2;

rNom=real(wi*a);  % nominal radius of shells, 2D function, as a function of theta
% Correction of membrne thickness for radius variations
rDeriv=real(wi*(1i*a.*(0:nterms-1)'));  % derivative of a wrt theta
thkFactor=sqrt((rDeriv./rNom).^2+1);

if rBoundFactor>0
    % Compute the smoothing function for angular dependences
    rBound=(rNom-mctr*thkFactor)*rBoundFactor;
    argSm=min(1,r./rBound);
    funSm=0.5*(1-cos(pi*argSm));
    wSm=w.*funSm;
    wiSm=wi;
    % wiSm(:,2)=wi(:,2).*funSm;
    for j=2:nterms
        wiSm(:,j)=wiSm(:,j-1).*wSm(:);
    end;
    rNomSm=real(wiSm*a);
else
    wiSm=wi;
    rNomSm=rNom;
end;

W1=zeros(N2,1,'single','gpuArray');
%     tan(alpha)=aDeriv; we want mbn thickness = 1/cos(alpha)
%     aDeriv=sin/cos so 1/cos=sqrt(aDeriv^2+1)
% rD2d=reshape(rDeriv,N1);
% rNom2d=reshape(rNom,N1);
for i=2:nm
    r0=rNomSm+(i-mctr-.5)*thkFactor;  % corrected radius as a fcn of theta.
    r02=r0.^2;  %max(0,reshape(r0,N1)).^2;
    b=model(i-1)-model(i);
    if abs(b)>tol*maxDensity  % A significant change in density
        srp2=sqrt(complex(r02-rp2));
        srm2=sqrt(complex(r02-rm2));
        W1=W1+b*real((rp.*srp2)+r02.*atan(rp./srp2)...
                     -rm.*srm2-r02.*atan(rm./srm2));
    end;
end;
W1=reshape(W1.*real(wiSm*aw),N1);  % impose the weighting

if any(N1<n)
    W=ExtractImage(W1,round(org),n,1);  % insert the image into the larger one.
else
    W=W1;
end;
