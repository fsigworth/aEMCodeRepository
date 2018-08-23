function [mf c T mfg]=meFilterImageFlatLF(mc,mi,LFAmp)
% function [mf c T mfg]=meFilterImageFlatLF(mc,mi,LFAmp)
% Boost low frequencies of the merged image mc to give it an effective ctf
% that is flat between f=0 and the first peak of the low-defocus ctf.  The
% image therefore approximates a phase-plate image.  c is the starting
% effective ctf (2D) and T is the transfer function of the filter.
% Called without arguments, this function puts up a file selector to load
% mi, and uses information in mi to load mc.  The output mfg is a
% gaussian-filtered version of mf.

if nargin<3
    LFAmp=0.4;
end;
if nargin<2  % Not enough arguments, go get the files.
    [infoname pa]=uigetfile('*.mat','Select an mi file');
    if ~ischar(pa)
        return
    end;
    
    cd(pa);
    
    load(infoname);                        % Get the mi structure
    imageName=[mi.basePath mi.procPath mi.baseFilename 'm.mrc'];
    mc=ReadMRC(imageName);  % Get the merged file
else
    imageName='';
end;

fc=.3;  % some Gauss filtering, 60% of Nyquist

%%
% figure(1); SetGrayscale;
% imacs(GaussFilt(mc,fc));
% title(imageName);
% drawnow;

n=size(mc);
ds=mi.imageSize(1)/n(1);  % current pixel size
pixA=mi.pixA*ds;
df=1/(pixA*min(n));       % frequency step, 1 pixel in short direction
c=meGetEffectiveCTF(mi,n); % assume we are filtering an entire image.
% Compute where the farthest first ctf zero could be.  The second zero will
% be 1.4 x farther out.
% minDefocus=mi.ctf(1).defocus-mi.ctf(1).deltadef;
minDefocus=mi.ctf(1).defocus;  % ignore astigmatism
alpha=mi.ctf(1).alpha;
lambda=mi.ctf(1).lambda;
nf0=sqrt(abs((1-alpha)*1e-4/(minDefocus*lambda)))*pixA; % norm. radius of 1st zero
nf1=sqrt(abs((0.5-alpha)*1e-4/(minDefocus*lambda)))*pixA; % norm. radius of 1st maximum
r=RadiusNorm(n);  % normalized frequency
mx1=mean(c(abs(r-nf1)<df));  % mean value of peak, averaged around a ring
dr=nf0*n(1)/5;  % feathering;
c1=c+(1-fuzzymask(n,2,nf0*n,dr));  % c out to first zero, then c+1 elsewhere.
T=(1-fuzzymask(n,2,nf0*n,dr))...   % 1 beyond the first zero
    + LFAmp/mx1*(fuzzymask(n,2,nf0*n,dr)-fuzzymask(n,2,nf1*n,dr))... %const in annulus
    + LFAmp*fuzzymask(n,2,nf1*n,dr)./c1;  % low frequency comp
% T=(1-fuzzymask(n,2,nf1,dr))...   % 1 beyond the first peak
%     + LFAmp*fuzzymask(n,2,nf1,dr)./c1;  % low frequency comp
% T= LFAmp/mx1*(fuzzymask(n,2,nf0,dr)-fuzzymask(n,2,nf1,dr))... %const in annulus
%     + LFAmp*fuzzymask(n,2,nf1,dr)./c1;  % low frequency comp

mf=real(ifftn(fftn(mc).*fftshift(T)));
%%
outName=[mi.basePath mi.procPath mi.baseFilename 'mlf.jpg'];

% figure(1);
% imacs(mfg);
% title(outName);

% figure(2);
% x=(0:n/2-1)*df;
% plot(x,[sectr(c) sectr(T.*c)]);
% xlabel('Frequency, A^{-1}');
% ylabel('CTF');

if nargin<2
imwrite(uint8(imscale(mf)),outName);
    mf=[];
end;

%%
mfg=GaussFilt(mf,fc);

figure(1);
imacs(GaussFilt(mfg,.15));