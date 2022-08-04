% FTEx5.m  Demonstrate power spectra
% Must have an image ms to read.
% load ms -mat
baseName='sq10_001';
ms=ReadEMFile([baseName 'ala.mrc']);
ms(:,:,2)=ReadEMFile([baseName 'alb.mrc']);
mi=ReadMiFile([baseName 'mi.txt']);
[n n1 n2]=size(ms);
df=1/mi.pixA/n;
figure(2); clf;
set(gcf,'color',[.4 .4 .4]);

% SetGrayscale;
nd=n/4;
bf=4;  % radial spectrum bin factor
dexp=.02;
ds=2;  % downsampling for radial spectrum
fs=[-nd/2 nd/2-1]*df;
defoci=[1.6 7.9];

for i=1:2
    subplot(2,3,3*i-2);  % show the filtered images
    m=ms(:,:,i);
    mfilt=GaussFilt(m,0.05);

    imags(mfilt); axis off;
        title([num2str(defoci(i)) ' um']);
        drawnow;

%     subplot(231);
%     m=ms(:,:,i);
%     imags(m); axis off; drawnow;
%     title([num2str(defoci(i)) ' um']);
    m=m-mean(mean(m));
    mf=fftn(m);
    sp=Crop(abs(fftshift(mf)),nd);
%     sp=sp(n/4+1:3*n/4,n/4+1:3*n/4);
    subplot(2,3,3*i-1);
    imags(fs,fs,(abs(GaussFilt(sp,0.05)).^dexp)); drawnow;
    rsp=RadialPowerSpectrum(BinImage(m,ds),0,4);
    subplot(2,3,3*i);
    semilogy(df*bf*(1:n/(2*bf*ds)),rsp);
    axis([0 .15 3 1000]);
    grid on;
    xlabel('Spatial frequency, A^{-1}');
end;


    