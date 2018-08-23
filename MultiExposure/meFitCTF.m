function [ctf, spectrum, disData]=meFitCTF(m,mi,pars,ctfOptions,imgIndex,doRemoveSpots)
% function [ctf spectrum, disData]=meFitCTF(m,mi,pars,ctfOptions,imgIndex,doRemoveSpots)
% Given an image m and parameter structures, returns the ctf
% structure.
% Also return the averaged spectrum in the struct:
% spectrum.df        % frequency step in A^-1
% spectrum.spectrum  % square 2D spectrum, f=0 in center
% defocus is the approximate defocus value.

if nargin < 6
    doRemoveSpots=0;
end;

% Get the image size
[nx, ny, nim]=size(m);
n=[nx ny];

% downsample if appropriate
maxRes=ctfOptions.maxRes(imgIndex);
minRes=ctfOptions.minRes(imgIndex);

ds=mi.imageSize(1)/nx;
while ds*4*mi.pixA<maxRes% at most half-nyquist
    ds=ds*2;
end;

ns=mi.imageSize/ds;
mc=Downsample(m,ns);  % reduce image size by ds
% pixA=mi.pixA*ds;

% CT fitting parameters
Pa.lambda=EWavelength(mi.kV);
Pa.theta=pi/180*(0:10:80);
Pa.alpha=ctfOptions.alpha;
Pa.Cs=2;  % 2 mm
Pa.B=ctfOptions.B;
% Pa.res=pixA;
defs=pars.searchDefoci(:,imgIndex);
Pa.defocus=defs(1):ctfOptions.defocusStep(imgIndex):defs(2);
dfrange = max(sqrt(Pa.defocus(1))/2,sqrt(Pa.defocus(end)/10));
Pa.deltadef=-dfrange:dfrange/6:dfrange;

% if doRemoveSpots
%     %         ns1=n/ds1;
%     nx=max(ns);
%     if nx>min(ns)  % rectangular image
%         win=SquareWindow(ns);
%         mcpad=Crop(win.*mc,nx);  % downsampled, padded image
%         fmc1=fftn(mcpad);
%     else
%         fmc1=fmc;
%     end;
%     fmcs=fftshift(fmc1);  % zero-centered
%     ds1=max(ds,2);
%     nr=nx/ds1;
%     if ds1>ds
%         fmcs=Crop(fmcs,nr);
%     end;
%     minr=nr*pixA*ds1/70;
%     threshSD=5;
%     display=pars.useGraphics;
%     %         disp('RemoveSpots');
%     [spc, pmask]=RemoveSpots(abs(fmcs).^2,minr,threshSD,display);
%     % pmask is now zeros at the spots, and square.
%     xpmask=Crop(pmask,nx,0,1);  % fill the outside with ones
%     rmask=Downsample(xpmask,ns);  % Convert to rectangle
%     mc=real(ifftn(ifftshift(rmask).*fmc));
% end;
%     disp(['CTF fitting ' num2str(i)]);

[P, c, sp, disData]=FitCTF(mc,Pa,mi.pixA*ds,maxRes,minRes,ctfOptions.showGraphics,ctfOptions);
drawnow;

spectrum.spectrum=sp;
spectrum.df=1/(size(sp,1)*mi.pixA*ds);
% % P.res=mi.pixA;  % Change it back to the original pixel size
% P.res=[];  % Change it back to the original pixel size
% P=rmfield(P,'pixA');
P.ampFactor=1;
ctf=P;


