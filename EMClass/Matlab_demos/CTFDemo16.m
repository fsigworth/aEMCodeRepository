% CTFDemo2
%  attempt at 
n=1024;
lambda=16;
a=1;
a0=1;
ctr=[0 0];
cz=n/2+1;

% Basic waves
[X,Z]=ndgrid(-n/2:n/2-1);

m=exp(1i*2*pi*Z/lambda);

R=sqrt(X.^2+Z.^2);
Rn=R;
Rn(R<1)=1;
s=1i*exp(1i*2*pi*R/lambda)./Rn;

figure(1);
SetComplex;

cz=n/2+1;

%%
for i=0:200
% imags(X(:,1),Y(1,:),abs(m+s*a).^2);
z=cz+i;
% figure(1);
% psi=a*s+a0*m;
% A=real(conj(m).*psi);
% imacx(a*s+a0*m);

figure(2);
pause(0.1);
I=squeeze(sqrt(abs(a*s+a0*m).^2)-a0);
% imags(A);
plot((A(:,cz+i)-1).*R(:,cz));
plot((A(:,cz+i)-1));
drawnow;
end;

