function out=BFilter(in, B, fc,fw)


if nargin<4
    fw=fc/10;
end;

b=B/4;
k=1.782/fw;

n=size(in,1);
[x,y,z]=ndgrid(-n/2:n/2-1);
r2=ifftshift(x.^2+y.^2+z.^2)/n^2;
r=sqrt(r2);  % normalized frequency

h=0.5*exp(-b*r2).*(1-erf(k*(r-fc)));
out=real(ifftn(h.*fftn(in)));

