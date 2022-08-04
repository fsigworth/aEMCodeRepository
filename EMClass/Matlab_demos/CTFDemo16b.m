% CTFDemo2
%  attempt at line diffraction
n=1024;
ny=128;
sz=[n ny n];

lambda=8;
a=.1;
a0=1;
ctr=sz/2+1;
cx=n/2+1;
cy=ny/2+1;

% Basic waves
[X,Y,Z]=ndgrid(-n/2:n/2-1,-ny/2:ny/2-1,1:n);

m=exp(1i*2*pi*Z/lambda);

R=sqrt(X.^2+Y.^2+Z.^2);
Rn=R;
Rn(R<1)=1;
s=1i*exp(1i*2*pi*R/lambda)./Rn;

figure(1);
SetComplex;


%%
i=0;
% for i=0:200
% imags(X(:,1),Y(1,:),abs(m+s*a).^2);
z=cx+i;
% figure(1);
psi=a*s+a0*m;
A=real(conj(m).*psi);
Af=GaussFilt(A,.01);
% imacx(a*s+a0*m);

% figure(2);
% pause(0.1);
% I=squeeze(sqrt(abs(a*s+a0*m).^2)-a0);
% If=GaussFilt(I,.02);
%%
imags(Af(:,cy,:));
% plot((A(:,cx,cx+i)-1));
% axis([0 inf -.5 .5]);
drawnow;
% end;

