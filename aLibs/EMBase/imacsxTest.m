% imacsxTest.m
% Test the new imacsx function

[mr, mi]=ndgrid(-.5:1/1024:511/1024);
imags(mr)
imags(mi)
mc=mr+1i*mi;
imacsx(exp(1i*mi*100))

cdddphi=0*mr;
amp=1;
phi(mi<0)=amp*cos(mr(mi<0)*200);
imacsx(exp(1i*(mi*1000+phi)))

