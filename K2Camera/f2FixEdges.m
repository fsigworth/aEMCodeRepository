function out=f2FixEdges(in)
% function out=f2EdgeFixer(in)
% replaces edge artifacts from a Falcon Hack image with constant pixels
xPix=69:3989;
yPix=98:4010;
n=size(in);
out=ones(n)*median(in(:));
out(xPix,yPix)=in(xPix,yPix);

