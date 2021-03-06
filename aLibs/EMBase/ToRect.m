 function W = ToRect(P, n, rscale, center)
%  function W = ToRect(P, n, rscale, center)
% W = ToRect(P, n, rscale, center)
% Convert the polar representation P(r,theta) to rectangular coords
% W, with the origin given by the vector center (default: n/2+1,n/2+1).
% The optional rscale parameter (default is 1) is the same as used in ToPolar; it is
% the scaling of r values to distance in W.  A value of 0.5 means
% that one unit in the r direction of P corresponds to an increment of 0.5
% in the distance from the center of W.
% If size(P,1)=1, we get only angular variations. if size(P,2)=1, we get
% only radial variations (rot symmetrized).

% % test code
% setgrayscale;
% figure(1);
% n=64;
% a=18;
% d=4;
% x=32;
% y=32;
% W=vidensity(n,a,d,x,y)+randn(64);
% % W(32:34,:)=-5;
% W(:,32)=-5;
% W(40,40)=-5;
% subplot(2,2,1);
% imagesc(W);
% axis xy;
%  
% ntheta=192;
% nr=32;
% rscale=1;
% x0=28;
% y0=28;
% P = ToPolar(W, ntheta, nr, rscale, x0, y0);
% subplot(2,2,2);
% imagesc(P);
% axis xy;
% n=64;
% center=[x0 y0];
% % end of first test code.

[nr, nt]=size(P);
if nr==1 % have to have at least two radial values.
    nr=2;               % (but one angular value is ok.)
    P=repmat(P,2,1);
end;

if nargin <3
    rscale=1;
end;
if nargin < 2
    n=2*nr;
end;
if nargin<4
    center=[n/2+1 n/2+1];
end;

eps=0.0001;
n2=n.^2;

% Expand the polar representation by one unit in theta, wrapping it.
P(:,nt+1)=P(:,1);

% Find the polar coordinates corresponding to each rectangular matrix element.
[X, Y]=ndgrid(1-center(1):n-center(1), 1-center(2):n-center(2));
R=1+sqrt(X.^2+Y.^2)/rscale;
T=nt/(2*pi)*atan2(Y,X)+1;
T=T+(T<1)*nt;  % Unwrap theta values.

% Convert to 1D arrays
R=reshape(R,n2,1);
T=reshape(T,n2,1);

% Clip the coordinates to be in bounds.
R=max(R,1+eps);
R=min(R,nr-eps);

T=max(T,1+eps);
T=min(T,nt+1-eps); % This value should never be exceeded anyway...

R0=floor(R);
Ri=R-R0;
R1=R0+1;

T0=floor(T);
Ti=T-T0;
T1=T0+1;
T0=(T0-1)*nr;  % Convert to linear index component for columns.
T1=(T1-1)*nr;

W= (1-Ti).*(P(R0+T0).*(1-Ri) + P(R1+T0).*Ri)...
     + Ti.*(P(R0+T1).*(1-Ri) + P(R1+T1).*Ri);
W=reshape(W,n,n);

% %  ----more test code follows----
% subplot(2,2,3);
% imagesc(W);
% axis xy;
