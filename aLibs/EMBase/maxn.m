function [val,inds]=maxn(m)
% function [val,inds]=maxn(m)
% find the maximum value and its indices in an n-dimensional array.
%
s=size(m);
[val,i]=max(m(:));
[a b c d e f g h]=ind2sub(s,i);
inds=[a b c d e f g h];
inds=inds(1:ndims(m));