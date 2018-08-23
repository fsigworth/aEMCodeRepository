% CT fitting parameters

mc=randn(3840);
pixA=1.78;

mi.kV=300;
ctfOptions=struct;
ctfOptions.alpha=.02;
ctfOptions.B=30;
defs=.5 :.2 :6;
% dfrange=.5;
showGraphics=1;

Pa.lambda=EWavelength(mi.kV);
Pa.theta=pi/180*(0:10:80);
Pa.alpha=.02;
Pa.Cs=2;  % 2 mm
Pa.B=ctfOptions.B;
% Pa.res=pixA;
% defs=pars.searchDefoci(:,imgIndex);
Pa.defocus=defs;
dfrange = max(sqrt(Pa.defocus(1))/2,sqrt(Pa.defocus(end)/10));
Pa.deltadef=-dfrange:dfrange/6:dfrange;

maxRes=5;
minRes=15;

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

[P, c, sp, disData]=FitCTF(mc,Pa,pixA,maxRes,minRes,showGraphics,ctfOptions);
drawnow;
