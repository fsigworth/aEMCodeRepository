% krLiveCTF.m

d=dir;
name=d(3).name
m=ReadMovie(name);
ms=RemoveOutliers(sum(single(m),3));
n=NextNiceNumber(size(m,2))
md=Downsample(Crop(ms,n,0,mean(ms(:))),n/2);

%%
figure(1);
mysubplot(111);
imaga(imscale(Downsample(md,n/4),256,.0001));
axis equal off;
title(name,'interpreter','none');
figure(2);
%%
opts=struct;
  Pa.lambda=EWavelength(300); % lambda in angstroms
  Pa.defocus=0:0.5:20;  % in microns
  Pa.deltadef=-.5:.1:.5;  % astigmatism
  Pa.theta=[-90:30:90]*pi/180;
  Pa.alpha=0:pi/10:pi/2;  % This can also be fitted, for use with phase plate.
  Pa.B=40;  % This is fitted if a vector of values is given.
   % The following parameters are not fitted
  Pa.Cs=2;
  P=FitCTFk2(md,Pa,1.05,2,[.05 .25],opts);

