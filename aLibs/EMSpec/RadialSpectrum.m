function [s,freqs]=RadialSpectrum(in,ctPars,pixA,binFactor)
% [s,freqs]=RadialSpectrum(in,ctPars,pixA,binFactor)
% Create the 1D average spectrum using astig and Cs corrections
% Based on RadialPowerSpectrum

if nargin<3 && isfield(ctPars,'pixA')
    pixA=ctPars.pixA;
end;
if nargin<4
    binFactor=1;
end;
n=size(in);
% Make a window of width n/8 around the edges.
w=SquareWindow(n, ceil(n(1)/8));
norm=w(:)'*w(:);  % Accounts for power attenuation due to the window.
ws=sum(w(:));

% Remove DC component.
in=in-mean(in(:));

in=in.*w;  % Windowed signal
xs=sum(in(:));
in=in-w*(xs/ws); % Second round of DC subtraction.

% Spectrum in rectangular coordinates
fs=fftshift(abs(fftn(in)).^2)/prod(n);
if binFactor>1
    fs=BinImage(fs,binFactor);
end;
% Radial averaging
[s,freqs]=RadialCTF(fs,ctPars,pixA);

