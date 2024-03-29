function m = fuzzymask(n,dims,r0,risetime,origin)
% function m = fuzzymask(n,dims,r0,risetime[,origin]);
% Create a centered 1d, 2d or 3d disc of radius r0, made with an error
% function with effective risetime tr (feathering width) in pixels.  By
% default risetime = r0 / 10.  m is double unless either n or r0 are
% singles.
% 
% If risetime=0 then the m is logical ones in all pixels with r<r0.  r0 is then the
% bounding outer radius of the disc.  To have a
% maximal disc with max raidus 10, fuzzymask(21,2,11,0); a disc with only
% single pixels at the edges is gotten with fuzzymask(21,2,10.001,0).
% 
% The old fuzzydisc has an fc argument that corresponds to tr in this way:
% fc=0.1 <--> tr=3.3
% fc=0.2 <--> tr=1.7
% The origin variable is a vector; its default value is floor(n/2)+1 in each
% dimension.
% - If dims=2 and r0 is a vector, an elliptical disc is made.  r0(1) is the
% radius along x; r0(2) is radius in the y direction.  r0=[rr rr]
% creates a circle; r0=[rr 2*rr] creates an ellipse with the major axis
% along y with radius rr*2, minor axis along x with radius rr.
% - If dims=3 and r0 is a vector, then an ellipsoid is made.  r0(1) is the
% radius rx along x; r0(2) is ry and r0(3) is rz.
% Note that the fuzzy boundary gets extended in y (2D and 3D case) and z
% (3D case) proportionally to the radii.

ctr=floor(n/2)+1;  % default center is appropriate for ffts.
if nargin<4
    risetime=r0/10;
end;
if risetime(1)==0
    k=0;
else
    k=1.782./risetime;
end;

switch dims
    case 1
        if nargin<5
            origin = ctr;  % Coordinate of the center
        end;
        r=abs(1-origin(1):n-origin(1))';

    case 2
        if numel(n)<2
            n(2)=n(1);
        end;
        if nargin<5
            origin = floor(n/2)+1; % Coordinates of the center
        end;
        [x,y]=ndgrid(1-origin(1):n(1)-origin(1),1-origin(2):n(2)-origin(2));
        if numel(r0)<2
            r = sqrt(x.^2 + y.^2);
        else
            r=sqrt(x.^2+(y*r0(1)/r0(2)).^2);
        end;

    case 3
        if nargin<5
            origin = [ctr ctr ctr]';  % Coordinate of the center
        end;
        [x,y,z]=ndgrid(1-origin(1):n-origin(1),1-origin(2):n-origin(2), 1-origin(3):n-origin(3));
        if numel(r0)<3
            r = sqrt(x.^2 + y.^2 + z.^2);
        else
            r = sqrt(x.^2 + (y*r0(1)/r0(2)).^2 + (z*r0(1)/r0(3)).^2);
        end;
    otherwise
        error(['Number of dimensions out of bounds: ' num2str(dims)]);
end;

if k==0
    m=r<r0;
else
    m=0.5*(1-erf(k(1)*(r-r0(1))));
end;