function out=Upsample2(in);
% function out=Upsample2(in);
% Increase the size of a 2D map by a factor of 2 in each dimension, copying pixels.
% fs 5 Dec 03

[nx ny]=size(in);

out(1:2:2*nx,1:2:2*ny)=in;
out(2:2:2*nx,1:2:2*ny)=in;
out(2:2:2*nx,2:2:2*ny)=in;
out(1:2:2*nx,2:2:2*ny)=in;
