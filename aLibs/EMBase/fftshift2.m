function out=fftshift2(in)
% simple shifting of a stack of images.
sz=size(in);
sz(3:end)=0;  % operate only on 1st two dimensions
out = circshift(in,floor(sz/2));  % code from fftshift
