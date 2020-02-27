% imacsxTest.m
% Test the new imacsx function
n=1024;
n2=n/2;
fx=100; % cycles per width
fy=100; % cycles per height
amp=.5; % extent of grating phase shift, radians

gratingY=.95; % height of grating
gratingX=.05; % half width of grating in x
xs=(-.5:1/n:.5-1/n)';
[mr, mi]=ndgrid(xs,-gratingY:1/n:1-gratingY-1/n);
imags(mr)
imags(mi)
mc=mr+1i*mi;
psi0=exp(1i*mi*2*pi*fy);
% imacsx(psi);
n0=round(gratingY*n)+1;  % zero-phase point

p=abs(xs)<gratingX;
phi=zeros(n,1);
phi(p)=amp*cos(xs(p)*2*pi*fx);
wave=exp(1i*phi);
% plot(phi)
% psi(:,n0)=phi;
Psi=zeros(n,n);

Psi(:,n0)=fft(wave);
dfz=.1;
for i=1:n0-1
    Psi(:,n0-i)=Psi(:,n0-i+1).*exp(-1i*(pi*dfz*(xs.^2)+2*pi*fy/(n+10)));
%     psi(:,n0+i)=ifft(Psi(:,n0+1));
end;
psi=ifft(Psi,[],1);
psi(:,n0+1:n)=psi0(:,n0+1:n);
mysubplot(121);
imacsx(psi-psi0);
mysubplot(122);
imags(abs(psi))


% phi=0*mr;
% amp=1;
% 
% p=mi<0 & abs(mr)<gratingX;
% phi(p)=amp*cos(mr(p)*200);
% imacsx(exp(1i*(mi*1000+phi)))

