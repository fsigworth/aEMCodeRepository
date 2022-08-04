% CTFDemo16c
%  wave illustration

d=5;
nd=5;
zLMin=-10;
zLMax=40;
w=d*nd;
theta=asin(1/d);
c=cos(theta);
s=1/d;
R=[c s; -s c];


figure(1)
xs=[-w w];
plot(xs,[0 0],'k-','linewidth',3);
hold on;
xds=[-nd*d:d:nd*d];
plot(xds,0*xds,'r.','markersize',30);

xs1=.1*xs;
lx=d*1.5;
% lx=d*3;
line=[-lx lx ; 0 0];
shft=[lx lx ; 0 0]*1;
% shft=[lx lx ; 0 0]*0;
cc0=.4;
cc=.4;  % color constant
cc2=mean([1 cc]);
colors=[cc0 cc0 cc0; 1 cc 1; 1 cc2 cc];
for zL=zLMin:zLMax
    p0=line+[0 0; zL zL];
%     p90=line+[0 0; zL+.25 zL+.25];
    p90=line+[0 0; zL+.75 zL+.75];  % correct: phase is advanced
    
    p1=R*(p90-shft);
    p2=inv(R)*(p90+shft);
    plot(p0(1,:),p0(2,:),'-','color',colors(1,:));
       plot(p1(1,:),p1(2,:),'-','color',colors(2,:));
    plot(p2(1,:),p2(2,:),'-','color',colors(3,:));
end;
zMaxFrac=.8;
lines=[0 0; 0 zLMax*zMaxFrac];
lines0=[0 0; -(s)*shft(1,1) (zLMax*zMaxFrac-s*shft(1,1))/c];
% lines=[0 0; zLMin*.1 zLMax*.3];
lines(:,:,2)=(R)*(lines0-shft);
lines(:,:,3)=inv(R)*(lines0+shft);
widths=[3 2 2];
for i=1:3
    plot(lines(1,:,i),lines(2,:,i),'-','color',colors(i,:),'linewidth',widths(i));
end;
plot(xs,[0 0],'k-','linewidth',3);
plot(xds,0*xds,'r.','markersize',20);
% plot(0,38,'ko','markersize',70,'color',[.5 .5 1]);
axis([-w w -10 50]);
% axis([-w/2 w/2 -5 25]);
hold off;

    



return

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

