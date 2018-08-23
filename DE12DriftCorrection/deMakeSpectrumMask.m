function [msk mskx msky]=deMakeSpectrumMask(sz, maskWidth, ds)
% function [msk mskx msky]=deMakeSpectrumMask(sz, maskWidth)
% Make a rectangular mask with relative thickness maskWidth around the
% edges of an image of size sz.  By default maskWidth = 0.1.  The thickness
% will be a multiple of ds in width, by default 8.
% mskx includes left and right edges; msky are top and bottom strips
% between the edges.  msk is the union of the two (nonoverlapping) regions.

if nargin < 2
    maskWidth=0.1;    % Normalized edge thickness for mask.
end;
if nargin < 3
    ds=8;
end;
if numel(sz)<2
    sz=[sz sz];
end;
nx=sz(1);
ny=sz(2);

ds2=2*ds;

% Make a mask along left and right edges, mainly shows white part of noise
nax=ds2*round(nx*maskWidth/ds2);
mskx=zeros(sz);
mskx(1:nax,:)=1;
mskx(nx-nax-1:nx,:)=1;

% Mask along the top and bottom, used for Gaussian component
nay=ds2*round(ny*maskWidth/ds2);
msky=zeros(sz);
msky(nax+1:nx-nax-2,1:nay)=1;
msky(nax+1:nx-nax-2,ny-nay-1:ny)=1;

% entire region used for fit
msk=mskx+msky;
