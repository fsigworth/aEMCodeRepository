% EstCarbonSpectrum

basePath='/EMWork/Hideki/130629carbon/Micrograph/';
preName='2nd_00';
postName='.mrc';
nim=6;
iCamera=1;

modelSpectrum=CCDModelSpectrum2D(iCamera);  % handle DE-12 or CCD

m=[];
figure(1);
clf;
SetGrayscale;
for i=1:nim
    fullName=[basePath preName num2str(i) postName]
    [m0, pixA]=ReadEMFile(fullName);
    m(:,:,i)=RemoveOutliers(m0);
    imacs(BinImage(m(:,:,i),4));
    title(fullName,'interpreter','none');
    drawnow;
end;
nx=size(m,1);
mw=mePreWhiten(m,modelSpectrum);
sps=fftshift(abs(fft2(mw).^2));
%%
spf=GaussFilt(sps,.1,1);
sp1=[];
for i=1:nim
    sp1(:,i)=Radial(spf(:,:,i));
    sMean=mean(sp1(800:1200,i));
    sp1(:,i)=sp1(:,i)-sMean;
end;
%%
% for i=1:nim
%     sp1(:,i)=RadialPowerSpectrum(mw(:,:,i));
%     sMean=mean(sp1(800:1200,i));
%     sp1(:,i)=sp1(:,i)-sMean;
% end;
% 

%%
freqs=(0:nx/2-1)/(nx*pixA);
plot(freqs,sp1(:,[1 2]));

dy=8e12;
axis([0 inf -0.1*dy dy]);

