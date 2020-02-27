% ColorWheel
% Make a colorwheel for complex number display

nx=640;
a=256;
[x y]=ndgrid(-nx/2:nx/2-1);
[r theta]=RadiusNorm(nx);
msk=fuzzymask(nx,2,a,5);
z=2*r.*exp(1i*theta).*msk;

figure(4);
imacx2(z);
axis equal off;
