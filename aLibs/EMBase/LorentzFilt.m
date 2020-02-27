function [out,H]=LorentzFilt(in, fc, stack, exponent)
% GaussFilt: 1D, 2D or 3D lorentzian filter.
% out = LorentzFilt( in, fc,stack,exponent)
% [out,H]=LorentzFilt(in,fc,stack,exponent) % returns the kernel H too.
% Filter the input matrix to give the half-power frequency fc
% (in units of the sampling frequency).  Thus a 1 kHz filter on data
% sampled at 10 kHz results if fc = 0.1.
% The filter of the form H=sqrt(1/1+f^exponent) so the power spectrum is Lorentzian
% but there is zero delay. The default exponent is 2.
% The output matrix is the same size as the input; the filter
% uses ffts so the boundaries are periodic. 1D inputs are assumed to be
% column vectors.
% If stack>0 then the last dimension of in is taken to be the number of
% images, which are processed separately.  These 'images' can be 1, 2 or
% 3D.

% F. Sigworth
if nargin<4
    exponent=2;
end;
if nargin<3
    stack=0;
end;
m=size(in);
ndim=max(1,sum(m>1)); % number of non-singleton dimensions
if stack
    ns=m(ndim);
    ndim=ndim-1;
else
    ns=1;  % number of stacked images
end;
if abs(fc)<1e-9
    out=in;
    return
end;
% The unit of frequency in the various dimensions will be
% 1/m(1), 1/m(2), etc.  We want the output to be 1/sqrt(2) at fc, i.e.
% The output kernel will be 1./sqrt(1+(x/(m*fc)).^2);

if ndim==1
    n=m(1);
    x=abs(-n/2:n/2-1)';
    H=sqrt(1./(1+(x/(m(1)*fc)).^exponent));	% Lorentzian kernel
    
elseif ndim==2
    [x,y]=ndgrid(-m(1)/2:m(1)/2-1, -m(2)/2:m(2)/2-1);
    H=sqrt(1./(1+(abs(x)/(m(1)*fc)).^exponent ...
        +(abs(y)/(m(2)*fc)).^exponent ));	% Lorentz kernel
elseif ndim==3
    [x,y,z]=ndgrid(-m(1)/2:m(1)/2-1, -m(2)/2:m(2)/2-1, -m(3)/2:m(3)/2-1);
    H=sqrt(1./(1+(abs(x)/(m(1)*fc)).^exponent ...
        + (abs(y)/(m(2)*fc)).^exponent + (abs(z)/(m(3)*fc)).^exponent ));	% Lorentz kernel
else
    error('LorentzFilt: input matrix has dimension > 3');
end;
Hs=fftshift(H);

if ns==1
    fq=Hs.*fftn(in);
    if isreal(in)
        out=real(fftn(conj(fq)))/prod(m); % This saves memory, to use fftn instead.
    else
        out=fftn(conj(fq))/prod(m);
    end;
else
    out=zeros(m);
    if ndim==1
        for i=1:ns
            fq=Hs.*fft(in(:,i));
            if isreal(in)
                out(:,i)=real(ifft(fq));
            else
                out(:,i)=ifft(fq);
            end;
        end;
    elseif ndim==2
        for i=1:ns
            fq=Hs.*fftn(in(:,:,i));
            if isreal(in)
                out(:,:,i)=real(ifftn(fq));
            else
                out(:,:,i)=ifftn(fq);
            end;
        end;
    elseif ndim==3
        for i=1:ns
            fq=Hs.*fftn(in(:,:,:,i));
            if isreal(in)
                out(:,:,:,i)=real(ifftn(fq));
            else
                out(:,:,:,i)=ifftn(fq);
            end;
        end;
    end;
    
end;
