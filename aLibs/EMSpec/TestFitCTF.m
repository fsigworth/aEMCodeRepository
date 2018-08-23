% TestFitCTF.m


  pixA=2.6;
  n=1980;
  sigmaC=1;
  sigmaShot=1;
  
  Pa.lambda=EWavelength(300); % lambda in angstroms
  Pa.defocus=4;  % in microns
  Pa.deltadef=.4;  % astigmatism
  Pa.theta=pi/180*-50;           % radians
  Pa.alpha=.02;  % This can also be fitted, for use with phase plate.
  Pa.Cs=2;
  Pa.B=80;  % This is not fitted
  Pa.pixA=pixA;
  
  m0=randn(n,'single');
  mShot=randn(n,'single');
  f=RadiusNorm(n)/pixA;  % spatial frequencies
  sp0=sigmaC*(f.^2+(.1)^2)./(f.^2+(.005)^2);
  c0=CTF(n,Pa);
  c1=c0.*sqrt(sp0);
  m=real(ifftn(fftn(m0).*ifftshift(c1)))+sigmaShot*mShot;
  %%
  Pf=Pa;
  Pf.defocus=1:.2:6;
  Pf.deltadef=0:.05:.4;
  Pf.theta=-pi/2:.1:pi/2;
%   Pf.theta=0;
  P=FitCTF(m,Pf,pixA,5,25);
  P.thetaDeg=P.theta*180/pi;
P
%%
opts=struct;
opts.maxRes=5;
opts.minRes=25;
opts.kV=300;
opts.alpha=.02;
opts.B=40;
P3=meFitCTF(m,pixA,2,0,opts);
  P3.thetaDeg=P3.theta*180/pi;

P3
