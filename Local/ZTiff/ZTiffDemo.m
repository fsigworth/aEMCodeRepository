% ZTiffDemo.m
cd('/Users/fred/EMWork/Hideki/121210/Micrograph')
m1=ReadEMFile('sq02_10000.tif');
m10=RemoveOutliers(m1);
%%
% m10=randn(4096,4096,'single');
%%
figure(4);
semilogy(RadialPowerSpectrum(m10,0,4));
%%
figure(1);
pars=struct;
pars.snrRatio=100;
pars.lfCutoff=.07;
pars.displayOn=1;
pixA=3;
filename='sq02_10000z.tif';
[nTrunc,ts,fsModel,m2,rm]=WriteZTiff(m10,pixA,filename,pars);
%%
figure(2);
fc=.2;
w=2;
m10f=GaussFilt(m10,fc);
m10fl=m10f;
y0=3900;
m10fl(:,y0-w:y0+w)=max(m10f(:));
subplot(1,2,1);
imags(m10fl)
xlabel('X');
ylabel('Y');
subplot(1,2,2);
plot([m10(:,y0) rm(:,y0)-250 m10(:,y0)-rm(:,y0)]);
legend('Original','Compressed','Difference');
xlabel('X');
ylabel('Pixel value');
% subplot(3,2,4);
% plot([m2(:,y0) m2(:,y0)-round(m2(:,y0))]);

%%
figure(3);
nd=512;
wc=0;
subplot(1,2,1);
m2c=m2(1:nd,end-nd+1:end);
y0c=nd/2+1;
m2cl=m2c;
m2cl(:,y0c-wc:y0c+wc)=m2cl(:,y0c-wc:y0c+wc)+max(m2c(:))/4;
imags(m2cl)
subplot(1,2,2);
xs=1:nd;
plot(xs,round(m2c(:,y0c)),'.-',xs,m2c(:,y0c)-round(m2c(:,y0c)),'.');

return
%%
pars=struct;
pars.snrRatio=100;
pars.lfCutoff=.4;
pars.displayOn=1;
pixA=3;
filename='sq02_10000wz.tif';
[nTrunc,ts,fsModel,m2,rm]=WriteZTiff(m2,pixA,filename,pars);
