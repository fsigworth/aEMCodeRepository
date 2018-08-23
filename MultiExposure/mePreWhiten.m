function mw=mePreWhiten(m,modelSpectrum)
% function mw=mePreWhiten(m)
% Pre-whiten the stack of unbinned images m.  Typically the modelSpectrum is
% larger than 1 near zero frequency (1/DQE), but it is normalized here, so
% the low-frequency amplitude of the image is unchanged.

if nargin<2  % assume the US4000
    modelSpectrum=CCDModelSpectrum2D(1,4096);
end;
[nx, ny, nim]=size(m);
    modelSpectrum=Crop(modelSpectrum,[nx ny]);  % In case we are downsampling
% Needs normalizing:
modelSpectrum=modelSpectrum/modelSpectrum(nx/2+1,ny/2+1);  % Normalize it

h=fftshift(sqrt(1./single(modelSpectrum)));
    mw=single(zeros(nx,ny,nim));
for i=1:nim
    mw(:,:,i)=real(ifftn(fftn(m(:,:,i)).*h));
end;
