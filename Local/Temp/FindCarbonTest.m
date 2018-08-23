% FindCarbon.m
% load Temp/h.mat

f1=.02;
f2=.1;
bp1=SharpHP(SharpFilt(ri,f2),f1);
v1=GaussFilt(bp1.^2,.05);
subplot(222);
imacs(v1);
subplot(224);
plot(v1(:,200))
subplot(221);
imacs(h.ifImageComp);

fi=meCTFInverseFilter(ri,h.mi,1,0,0);

subplot(221);
imacs(fi);
subplot(222);
plot(sect(fi)-0*sect(sqrt(max(v1,0))));