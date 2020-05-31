% SpectrumCompare.m
% Compare 160909 with 191228 datasets

basePath='~/mini_drobo/';
basePath='/Volumes/Drobo4/';
path16=[basePath 'siggpu_data/160909/KvLipo121_2/'];
baseName1603='sq02_1_0014_Sep09_19.15.05'; % 3.02 um
baseName1602='sq02_1_0028_Sep09_19.26.19'; % 2.26um
info1602=ReadMiFile([path16 'Info/' baseName1602 'mi.txt']);
mic1602=ReadMRC([path16 'Micrograph_sq02_and_sq81/' baseName1602 'ala.mrc']);

%%
path19=[basePath '191228.1/'];
baseName1902='GridSquare_28970644_FoilHole_28972832_02'; % 2.26 um
info1902=ReadMiFile([path19 'Info/' baseName1902 'mi.txt']);
mic1902=ReadMRC([path19 'Merged/' baseName1902 'm.mrc']);

%%
mag=info1902.pixA/info1602.pixA;
% npix=NextNiceNumber(size(mic1902,1)*mag)
npix=2048;
msc19=DownsampleGeneral(mic1902,npix,mag);
msc16=Crop(mic1602,npix);
pixA=info1602.pixA;

figure(1);
mysubplot(221);
imags(msc16);
mysubplot(222);
imags(msc19);

sp16=RadialPowerSpectrum(msc16);
mysubplot(223);
semilogy(sp16);

sp19=RadialPowerSpectrum(msc19);
mysubplot(224);
%% Manual fit
% xpos are where n=f^2*lambda*def, i.e. zeros of the ctf
fmax=1/(2*pixA)
maxN=fmax.^2*22600*.025
lambda=EWavelength(200);
def=info1602.ctf(1).defocus*1e4;
xpos=sqrt((1:maxN)'/(lambda*def))*npix*pixA;

freq=Radius(npix);
freq1=(0:1023)';
a1=2;
a2=2.5;
a3=1.5;
a4=3;

sp19f=11000*sp19.*( (1-a1*1e-4*freq1).*(1+freq1.^2*a2*1e-7) ...
    .*(1+freq1.^3*a3*1e-10).*(1+freq1.^4*a4*1e-13) ).^2;
semilogy(min(200,[sp19f sp16]));
hold on;
plot(xpos,sp19f(round(xpos)),'k.','markersize',10);
hold off;


%% Simplex fit to minima of spectra

xpts=round(xpos);
zFreq=freq1(xpts);
% y19=p(5)*1e4*sp19(xpts);
y16=sp16(xpts);
scl=[-1e-4 1e-7 1e-10 1e-13];
p=[1 1 1 1 1];
p=Simplex('init',p);

for iter=1:1000
    h=1;
    for i=1:4
        h=h.*(1+p(i)*scl(i)*zFreq.^i);
    end;
    y19=h.^2.*sp19(xpts)*p(5)*1e4;
if mod(iter,10)==1
    semilogy(xpos,min(10000,[y19 y16]),'-.');
    drawnow;
end;
err=(log(y19)-log(y16))./(xpos+400);
    se=err'*err;
    p=Simplex(se);
    title(num2str(p));
end;
%% show the final result
    h=1;
    for i=1:4
        h=h.*(1+p(i)*scl(i)*freq.^i);
    end;
    h1=sectr(h);
    y19=h1.^2.*sp19*p(5)*1e4;
    semilogy(min(10000,[y19 sp16]),'-');
    drawnow;

%% Having matched the backgrounds, now filter the entire micrograph 19 to match 16

    
    
    