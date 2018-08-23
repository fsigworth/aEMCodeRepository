function [mc pmask]=meRemoveCrystal(m,pixA,threshSD,displayOn)
% function [mc pmask]=meRemoveCrystal(m,pixA,displayOn)
% Remove streptavidin crystal periodicity from image m.
% pixA is the pixel size in angstroms.  The optional argument threshSD is
% the threshold for masking of Fourier spots.  threshSD=6 makes
% good-looking images; the value 5 is better to prepare the image for CTF
% determination.
% This version should work with rectangular images, but hasn't been tested
% with them.  fs 20 Apr 12

if nargin<3
    threshSD=6;
end;
if nargin<4
    displayOn=0;
end;
% m=mc;
% pixA=2.7;
% threshSD=6;
% displayOn=1;

sz=size(m);
nx=max(sz);  % larger dimension
m=m-mean(m(:));
sqwin=SquareWindow(sz);
mpad=Crop(m.*sqwin,nx);  % image is padded to be square

ds=floor(5/pixA);  % We pick up data out to 10A, need sampling at least 5A.
ds=NextNiceNumber(ds,7,-2);  % Find the next lower
ds=max(ds,1);
ns=nx/ds;            % Downsampled image size
% create the cropped FT
fmp=fftn(mpad);
fmc=Crop(fftshift(fmp),ns); % crop the fft of the padded image.

if displayOn
    figure(1); clf;
    SetGrayscale;
    imacs(real(ifftn(fftshift(fmc))));
    drawnow;
end;

minr=nx*pixA/70;
% ThreshSD=6;
% disp('RemoveSpots');
[spc pmask]=RemoveSpots(abs(fmc).^2,minr,threshSD,displayOn);
xmask=1-Crop(1-pmask,nx);  % pad the mask to match square image
rmask=Downsample(xmask,sz); % change aspect ratio
mc=real(ifftn(fftshift(rmask).*fftn(m)));  % block the spots
if displayOn
    subplot(2,3,6);
    imacs(BinImage(mc,ds));
end;
