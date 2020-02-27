function filtMap=SolventAndProteinGridder(map,fc,df)
% function filtMap=SolventAndProteinGridder(map)
% Given map as an output of the SolventAndProteinDensity function, which
% has a Gaussian of sd=1 at each atom position, perform an inverse filter
% to provide a sharp cutoff of fc (default is 0.4).
if nargin<2
    fc=0.4;
end;
if nargin<3
    df=.1;
end;
n=size(map);
% Get the sharp filter kernel
[~,h]=SharpFilt(map,fc,df);
sigma=n/(2*pi);
hInv=h./ifftshift(Gaussian(n,3,sigma));
filtMap=real(ifftn(fftn(map).*hInv));
