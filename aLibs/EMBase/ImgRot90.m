function out=ImgRot90(in,k)
% Same as rot90, but handles even-sized images correctly with center at
% n/2+1.  Like rot90, rotates the first two dimensions of higher-
% dimensional arrays too.  So ImgRot90(v) rotates a volume v about the z axis
if nargin<2
    k=1;
end;
k=mod(k-1,4)+1;  % in the range 1..4
n=size(in,1);
even=(mod(n,2)==0);
shifts=[1 0; 1 1; 0 1; 0 0];

out=rot90(in,k);
if even  % even number of pixels
    out=circshift(out,shifts(k,:));
end;
