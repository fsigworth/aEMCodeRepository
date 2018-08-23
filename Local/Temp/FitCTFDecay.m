% FitCTFDecayB.m
% Experiment with K2 images of blank beam and thin carbon.
% 
%  cd('/Users/fred/EMWork/Hideki/140625/control/movie_frames/control')
% mRef=ReadEMFile('CountRef_Jun25_10.00.26.dm4');
% mv0=ReadMovie('Jun25_10.14.06_bright.tif');
% mv1=ReadMovie('Jun25_10.11.09_thin_carbon_u1000_OLA50.tif');
% %%
% m0=RemoveOutliers(flip(sum(single(mv0),3),2).*mRef);
% m1=RemoveOutliers(flip(sum(single(mv1),3),2).*mRef);
mode=3;
n=3840;

% Fitting spectral data from Li...Agard 2013
% doseRates=[1.023 2.246 4.0542 8.1623 9.7326 12.87 14.9 20.158];
% bVals=[.0048 .0432 .0991 .2223 .2714 .3631 .4226 .5701];
% fit to b as a funciton of doseRates:
% b=-.0274+.0315*doseRate-9.355e-5*doseRate.^2;
% spectrum = N(1 - b*sinc^2(f*2.584)); with N being total counts in image
% Detection fraction:
% detFractions=.866*exp(-doseRates*5.149/400);  % fraction of counts
% Model for spectrum
% f=(0:n/2-1)'/n;
% r=8;  % dose rate
% b=-.0274+.0315*r-9.355e-5*r^2;
% S=(1-b*sinc(f*2.58).^2);
% plot(f,S);
% return


switch mode
    case 1
cd('/Users/fred/EMWork/Hideki/140625/control/');
mi1=ReadMiFile('Info/control_001_Jun25_10.11.09_thin_carbon_u1000_OLA50mi.txt');
mi1=ReadMiFile('Info/control_001_Jun25_10.11.09_thin_carbon_u1000_OLA50mi.txt');
m1=ReadEMFile([mi1.imagePath mi1.imageFilenames{1}]);
mi0=ReadMiFile('Info/control_002_Jun25_10.14.06_brightmi.txt');
m0=ReadEMFile([mi0.imagePath mi0.imageFilenames{1}]);
    case 2
cd('/Users/fred/EMWork/Hideki/140618/Control/');
m1=RemoveOutliers(ReadEMFile('R_thin_carbon_u1000_8.25sec.mrc'));
m0=RemoveOutliers(ReadEMFile('R_bright_u1000_8.25sec.mrc'));
m0=Crop(m0,n,0,mean(m0(:)));
m1=Crop(m1,n,0,mean(m1(:)));
    case 3
% cd('/Volumes/WDTera2/EMWork/Hideki/140701/control');
cd('~/EMWork/Hideki/140701/control');
d=dir('Info/');
dnm=d(7).name
mi1=ReadMiFile(['Info/' dnm]);
m1=ReadEMFile([mi1.imagePath mi1.imageFilenames{1}]);
cnm=d(9).name
mi0=ReadMiFile(['Info/' cnm]);
m0=ReadEMFile([mi0.imagePath mi0.imageFilenames{1}]);
Pa.defocus=4.8:0.05:6.4;


end;

figure(1);
SetGrayscale;
imacs(GaussFilt(m0,.1));
n=size(m0,1);
%%
pixA=1.247;
sp0=RadialPowerSpectrum(m0);
sp1=RadialPowerSpectrum(m1);
%%
sp0(1:10)=sp0(10);
sp1(1:10)=sp1(10);

f=(1:numel(sp0))'/(2*numel(sp0)*pixA);
semilogy(f,[sp0 sp1]);
fo=.185;
fit0=85-20*exp(-(f/fo).^2.7);
plot(f,[sp0 fit0]);
% plot(f,[sp0 sp1 fit0]);
% plot(f,sp1-fit0);
% m1c=Crop(m1,3710);
ctr=n/2+1;
sp12=fftshift(abs(fftn(m1)).^2)/(n^2);
sp12(ctr-10:ctr+10,ctr-10:ctr+10)=fit0(10);
fit2=ToRect(fit0);
%%
q=GaussFilt(sp12-fit2,.05);
imacs(abs(q).^.3);

%%
%
sc=sp12-fit2+fit2(1);
mc=m1;
pixA=1.247;
maxres=4;
minres=12;
Pa.lambda=EWavelength(200);
% Pa.defocus=0.8:0.1:1.5;
Pa.deltadef=-.5:.1:.5;
Pa.theta=0;
Pa.alpha=0:.03:.12;
Pa.Cs=2;
Pa.B=0;
test=1;
% [P, c]=FitCTF(mc,Pa,pixA,maxres,minres,test);

P=mi1.ctf;
%%
f0=(0:n/2-1)'/n;
mtf=CCDModelMTF(f0,5);
P.B=55;

P1=P;
P1.B=P.B/4;
% P1.deltadef=0;  % astigmatism makes no change in 1D spectrum.
P1.alpha=.05;
B=P.B;
ct=CTF(n,pixA,P1);
subplot(122)
imacs(ct.^2);
subplot(121);
imacs(sc.^.2);
cr=mtf.^2.*Radial(ct.^2);
sr=Radial(sc);
f=(0:n/2-1)'/(n*pixA);
% k=.15;
% k4=.029;
% a4=.0009;
a0=50;
k4=.034;
a4=2080;
expo=3.5;
subplot(1,1,1);
env=(a4*k4^expo)./(f.^expo+k4.^expo)+a0;
bEnv=exp(-B/2*f.^2);  % B4 but squared
fitC=1*cr.*env+80;
plot(f,[sr fitC env*1.3 bEnv*400 mtf.^2*400]);
axis([0 inf 0 800]);
legend('Spectrum','Fit','Carbon model',['B factor=' num2str(B)],'MTF^2');
title(dnm,'interpreter','none');

