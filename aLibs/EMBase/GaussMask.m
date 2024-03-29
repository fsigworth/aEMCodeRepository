function m = GaussMask(n,dims,sigma,origin)
% function m = GaussMask(n,dims,sigma,origin)
% Create a centered 1d, 2d or 3d Gaussian of s.d. sigma.  If sigma is a
% scalar, m is a Gaussian disc or sphere; if sigma is a vector, it is a
% Gaussian ellipse/ellipsoid.  The function is non-normalized, m(origin)=1.
% For 1d, m is returned as a column vector.
% To make a Gaussian filter of corner frequency fc, given sampling
% frequency fs and no. of points n, set
% sigma = n*fc/fs/sqrt(log(2)).  Since 1/sqrt(log(2))=1.2011,
% sigma = 1.2011*n*fc/fs
% so the m-dimensional Gaussian filter function could be implemented as
% GaussFilt(data,fc)=real(ifftn(fftn(data).*GaussMask(n,m,1.2011*n*fc)))

if nargin<4  % Set the default center to be correct for even and odd FFTs.
    origin=floor((n/2+1))*ones(dims,1);
end;

if numel(sigma)==1  % If sigma is scalar, duplicate it.
    sigma=sigma*ones(dims,1);
elseif numel(sigma)~=dims
    error('GaussMask: dimension of sigma doesnt match dims');
end;

k=1./(sqrt(2)*sigma);

switch dims
    case 1
        r2=((k*(-origin(1)+1:n-origin(1))).^2)';

    case 2
        [x,y]=ndgrid(k(1)*(-origin(1)+1:n-origin(1)),k(2)*(-origin(2)+1:n-origin(2)));
        r2 = x.^2 + y.^2;

    case 3
        [x,y,z]=ndgrid(k(1)*(-origin(1)+1:n-origin(1)),...
            k(2)*(-origin(2)+1:n-origin(2)), k(3)*(-origin(3)+1:n-origin(3)));
        r2 = x.^2 + y.^2 + z.^2;
    otherwise
        error(['GaussMask: number of dims can only be 1..3, but is ' num2str(dims)]);
end;

m=exp(-r2);
