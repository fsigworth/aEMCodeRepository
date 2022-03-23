function m = disco(n, r)
% DISCO
% m = disco(n, r): create an n x n matrix that is all zeros except
%	for a offset circular region of radius r in the center, where elements
%	are equal to 1.  The origin is [n/2+1, n/2+1].

%	Construct the disc
n=n(1);
r0 = n/2+1;  % Coordinate of the center
ep = 1e-4;
m=zeros(n);
for i= ceil(r0-r+ep) : floor(r0+r-ep)	% Loop over rows
	y= i - r0;
	dx=sqrt(r*r-y*y)-ep;
	m(i, ceil(r0-dx):floor(r0+dx)) = 1;
end

% End disc.
	
