% PhaseFlipDemo

n=1024;
pixA=.5;
sigma=10; % Å
[X,Y]=ndgrid((-n/2:n/2-1)*pixA);
freqs=RadiusNorm(n)/pixA;
f0=.25;
msk=GaussMask(n,2,sigma/pixA);
msk=Crop(ones(5*sigma),n);
m=msk.*cos(2*pi*X*f0);
imags(m);

q=sect(m);
qs=cumsum(q.^2);
qs([1 400 600 1000])
qs(600)-qs(400)

d=1;
c=CTF(n,pixA,.025,d,2,0,.02);
c=repmat(sect(c),1,n);
c=abs(c);
% c=c.^2;
mf=real(ifftn(fftn(m).*ifftshift(c)));

imags(mf);
plot(sect(mf))
q=sect(mf);
qs=cumsum(q.^2);
plot(qs);

qs([1 400 600 1000])
qs(600)-qs(400)

% plot(sectr(freqs),sectr(c))