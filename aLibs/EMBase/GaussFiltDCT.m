function out=GaussFiltDCT(in, fc, isStack)
% function out=GaussFiltDCT(in, fc, isStack)
% GaussFiltDCT: 1D or 2D gaussian filter using the discrete cosine
% transform to avoid "wrap-around" edge effects.
% Filter the input matrix to give the half-power frequency fc
% (in units of the sampling frequency).  Thus a 1 kHz filter on data
% sampled at 10 kHz would have fc = 0.1.
% A typical value of fc for an anti-aliasing filter is 0.2.
% The output matrix is the same size as the input; the filter
% uses ffts so the boundaries are periodic.
% If stack>0 then the last dimension of in is taken to be the number of
% images, which are processed separately.  These 'images' can be 1, 2 or
% 3D.
% The corner frequency fc is related to the standard deviation sigma of the
% impulse response by the following:
%           fc=.133/sigma
% fc is related to a B factor (EM definition exp(-Bf^2)) according to
%  fc=sqrt(log(2)/(2*B))*pixA where pixA is number of angstroms per pixel,
%  and B has units of A^2.

% F. Sigworth
% DCT mode 1 Mar 12 derived from GaussFilt().

if nargin<3
    isStack=0;
end;
m=size(in);
ndim=sum(m>1); % number of non-singleton dimensions
if isStack
    ns=m(ndim);
    ndim=ndim-1;
else
    ns=1;  % number of stacked images
end;

% The unit of frequency in the various dimensions will be
% 1/m(1), 1/m(2), etc.  We want the output to be 1/sqrt(2) at fc, i.e.
% at fx=fc*m units.  The output will be exp(-(x.^2/fx.^2)*ln(2)/2);

out=zeros(m);
k=-log(2)./(2*double(fc)^2*m.^2);

switch ndim
    case 1       % there isn't a 1D dct in Matlab, so we make our own.
        n=m(1);
        x=(0:n-1)';
        k=k(1)/4;
        q=exp(k*x.^2);
        for i=1:ns
            fq=q.*dct1(double(in(:,i)));
            if isreal(in)
                out(:,i)=real(idct1(fq));
            else
                out(:,i)=idct1(fq);
            end;
             out(:,i)=(out(:,i)+[out(2:n); out(n)]+[out(3:n); out(n); out(n)])/3;  % rough correction for flaky dct
        end;
    case 2
        [x,y]=ndgrid(0:m(1)-1, 0:m(2)-1);
        q=exp(k(1)/4*x.^2+k(2)/4*y.^2);	% Gaussian kernel
        for i=1:ns
            fq=q.*dct2(in(:,:,i));
            %                out(:,:,i)=real(idct2(fq));
            %                 else
            out(:,:,i)=idct2(fq);
            %                 end;
        end;
    otherwise
        error('DCT mode supports only 1D and 2D data');
end;
