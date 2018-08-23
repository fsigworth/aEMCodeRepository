function d=EllipsoidDensity(n,b,e,org)
if nargin<4
    org=ceil((n(1)+1)/2)*[1 1];
end;
e
nume=numel(e);
e(nume+1:4)=0;
e
% n=256;
% b0=90;
% org=[127 127];
% e=.3;
nt=2*n;  % number of theta steps
er=e(1);
ei=e(2);
es=e(3);
et=e(4);

% Create the r field
[r thetas]=Radius(n,org);
itheta=round(1+(1+thetas/pi)*nt);
dtheta=pi/nt;



% a=b0./sqrt(1-real(e)*cos(2*theta)-imag(e)*sin(2*theta));




rm=r-.5;
rp=r+.5;

a2=b^2./(1-er*cos(2*dtheta*itheta)-ei*sin(2*dtheta*itheta)...
    -es*cos(3*dtheta*itheta)-et*sin(3*dtheta*itheta));
d=real(rp.*sqrt(a2-rp.^2)+a2.*atan(rp./(sqrt(a2-rp.^2)))...
    -rm.*sqrt(a2-rm.^2)-a2.*atan(rm./(sqrt(a2-rm.^2))))./sqrt(a2)*b;

% imacs(d./sqrt(a2)*b);

