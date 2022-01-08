function [S, cameraName]=CCDModelSpectrum2D(iCamera,n,ds)
% Return the analytical fit to the spectrum of a CCD camera
% as determined with MeasureCCDSpectrum2D.  iCamera is an index into a table
% of fitted values.  At present the indices are
% iCamera=1 : Yale F20 US4000
% iCamera=2 : Okazaki JEM 2200 TVIPS camera
% iCamera=3 : DE12
% iCamera=4 : BrandeisFalcon1
% iCamera=5 : K2
% iCamera=6 : RIKEN Falcon2
% (For all but iCamera=3, we just call CCDModelSpectrum with a 2D array of
% frequencies.  The SpectModel2D.mat data are exactly the same as the
% SpectModel.mat data except for iCamera=3.)
% n is a scalar or 2-element vector giving the size of the image to be
% returned. Default is the full camera size (e.g. 3k x 4k for DE12, but the
% padded size 3840x3840 for K2).
% If n is smaller than the camera size, then the image is assumed to be
% downsampled by ds = nativeSize(1)/n(1), and the edges represent f_Nyquist/ds = 0.5/ds.
% f=0 is in the center.

if nargin<1
    iCamera=1;  % default is CCD
end;
switch iCamera
    case 3
        n0=[3072 4096];
    case 5
        n0=[3840 3840];  % k2 camera
    case 7
        n0=[5760 4096]; % k3 camera
    otherwise
        n0=[4096 4096]; % indices 1,2,4
end;
if nargin<2
    n=n0;
end;
if nargin<3
    ds=n0(1)/n(1);  % possible downsampling
end;
if numel(n)==1 % square image
    n=[n n];
end;

% Find our local directory
pa=fileparts(which('CCDModelSpectrum2D'));
% Retrieve parameters from SpectModel2D.mat in the local directory
if numel(pa)<1
    pa='.';
end;
load([pa '/SpectModel2D.mat']);

if iCamera>numel(models) && iCamera~=7  % 7 is a special case
    error('CCDModelSpectrum2D: iCamera index too large');
end;
models(7)=models(5);  % special case

p=models(iCamera).spectPars;
cameraName=models(iCamera).camera;

switch iCamera
    case 3   % special case of the DE camera
        
        [f, theta]=RadiusNorm(n);
        f=f/ds;
        % Aliasing code copied from ccFit2DSpectrum
        ang=max(abs(cos(theta)),abs(sin(theta)));
        fAliased=1./ang-f;
        y0=ccDE12Model(f,theta,p);
        ya=ccDE12Model(fAliased,theta,p);
        S=y0+ya;
        
    case {5 7}  % K2 camera, assume a constant spectrum
        S=models(5).spectPars(4)*ones(n,'single');
        
%     case 6  % doesn't quite work yet.
%         [f, theta]=RadiusNorm(n);
%         f=f/ds;
%         % Aliasing code copied from ccFit2DSpectrum
%         ang=max(abs(cos(theta)),abs(sin(theta)));
%         fAliased=1.1./ang-f;
%         y0=CCDModelSpectrum(f,iCamera);
%         ya=CCDModelSpectrum(fAliased,iCamera);
%         S=y0+ya;
        
    otherwise  % we could do this for the K2 as well.
        f=RadiusNorm(n)/ds;
        S=CCDModelSpectrum(f,iCamera);
end;
