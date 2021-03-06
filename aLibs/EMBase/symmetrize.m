function out=symmetrize(in, symmetry);
% out = symmetrize( in, symmetry );
% Force n-fold symmetry by rotationally averaging the in matrix
% to form the out matrix.  Uses the mrotate convention of
% rotation about x,y = (n+1)/2.
%
out=in;
qrot = 2*pi/symmetry;
for i=2:symmetry
	out=out+mrotate(in,qrot*(i-1));
end;
out=out/symmetry;
