% FSCDemo.m

fscSum=zeros(96,1);
nSum=0;

for i=1:10
a0=.002;
s1=squeeze(m1(:,97,:));
s1f=s1+6*GaussHP(s1,.3)+a0*randn(192);
s2=squeeze(m2(:,97,:));
s2f=s2+6*GaussHP(s2,.3)+a0*randn(192);
f1=fftshift(fftn(ifftshift(s1)));
f1n=fftshift(fftn(ifftshift(s1f)));
f2n=fftshift(fftn(ifftshift(s2)));

figure(5);
SetComplex;
subplot(221);
imags(s1f);
axis equal off
title('Half 1');

subplot(222);
imags(s2f);
axis equal off
title('Half 2');

subplot(223);
imacx(Crop(f1n,128),.2);
axis equal off;

subplot(224);
imacx(Crop(f1n,128),.2);
axis equal off;

%%
corr=real(f1n.*conj(f2n))./(abs(f1n).*abs(f2n));
fs=(1:96)/(192*1.3);
figure(6);
plot(fs,Radial(corr))
drawnow;
fscSum=fscSum+Radial(corr);
nSum=nSum+1;

end;
fscSum=fscSum/nSum;

%%
plot(fs,fscSum,'b-');
hold on
plot(fs,fs*0,'k-');
plot(fs,fs*0+.143,'r--');
hold off;
xlabel('Spatial frequency, Å^{-1}');
ylabel('Correlation');