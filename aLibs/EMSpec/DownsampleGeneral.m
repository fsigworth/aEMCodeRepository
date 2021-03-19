function [out, finalmag, nix]=DownsampleGeneral(in,nout,mag,stack)
% function [out finalmag nix]=DownsampleGeneral(in,nout,mag,stack)
% Change magnification by an arbitrary ratio.  3D inputs are assumed to be
% volumes unless stack=1. If mag is not given, or given as [] or 0, it is
% set to nout(1)/size(in,1).  The function supports odd and even-sized
% square or cube inputs.
% 
% It works by first padding the input image to nix, then taking fft, then
% cropping this oversampled fft to nox, do ifft, crop to nout.  The
% magnification is therefore nox/nix.
% 
% maxOversample, which is by default 1.3 for large (>128 pixel) inputs and
% outputs, limits the maximum padding.  The returned value nix = size(in) x
% oversample. The finalmag value, the actual magnification ratio which
% approximates mag, is returned.  Errors are typically < 0.1%

ni=size(in,1);
ndim=sum(size(in)>1);

if nargin<3 || numel(mag)<1 || any(mag==0)
    mag=nout(1)/ni(1);
end;

% if nargin<4 || maxOversample<1
    nx=max(ni,nout/mag);
    if nx>128
        maxOversample=1.3;
    else
        maxOversample=256/(max(nx)+16);
    end;
% end;
if nargin<4
    stack=false;
end;

% Trivial case of no magnification change
if mag==1
    out=Crop(in,nout,stack);
    finalmag=1;
    nix=ni(1);
    return
end;

%

% Do a brute-force search
nos=nout(1):maxOversample*nout(1);  % vector of nos values to try
nis=round(nos/mag);        % optimal corresponding nis values
errs=(nos/mag-nis)./nis;   % relative errors
[~, ind]=min(abs(errs)); % find the index of the best nos and nis

% tol=9e-3;
% ind=find(abs(errs)<tol,1,'first');
% if numel(ind)<1
%     error('tolerance too small');
% end;

% subplot(212); plot(errs,'.-'); title(nis(ind));
% nix=nis(ind)
% nox=nos(ind)

nix=nis(ind);
nox=nos(ind);
finalmag=nox/nix;

sz=size(in);
szo=0*sz+nout;
if stack  % stack of 2d images
    if ndim>3
        error('Stacks of volumes not supported');
    end;
    nim=sz(3);
    szo(3)=nim;
    out=zeros(szo,'single');
    for i=1:nim
        inp=Crop(in(:,:,i),nix);  % pad the input
        fout=Crop(fftshift(fftn(ifftshift(inp))),nox);
        out(:,:,i)=Crop(fftshift(real(ifftn(ifftshift(fout)))),nout)*finalmag^ndim;
    end;
else
    inp=Crop(in,nix);  % pad the input
    fout=Crop(fftshift(fftn(ifftshift(inp))),nox);
    out=Crop(fftshift(real(ifftn(ifftshift(fout)))),nout)*finalmag^ndim;
end;
