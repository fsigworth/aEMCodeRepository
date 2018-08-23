% function out=InverseBFilter(in, pixA, B, fc, deltaF)
% InverseBFilter: 1D, 2D or 3D gaussian filter.
% out = InverseBFilter( in, B, fc, deltaF)
% Apply the inverse-Gaussian filter
% H=exp(Bf^2/4) where f is in units of A^-1

n=size(in);
nim=size(in,4);
n=n(1:3);

k=1./(fc^2*pixA^2*n.^2);  % sharp filter argument
b=B./(4*pixA^2.*n.^2);

        [x,y,z]=ndgrid(-floor(n(1)/2):ceil(n(1)/2-1),...
                       -floor(n(2)/2):ceil(n(2)/2-1),...
                       -floor(n(3)/2):ceil(n(3)/2-1));
        if deltaF<=0
            q=(k(1)*x.^2+k(2)*y.^2+k(3)*z.^2)<1; % totally sharp filter
        else
            r=sqrt(k(1)*x.^2+k(2)*y.^2+k(3)*z.^2);
            w=deltaF/fc;
            q=0.5*(1-erf((r-1)/w));  % error function of radius
        end;
        q=q.*exp(b(1)*x.^2+b(2)*y.^2+b(3)*z.^2);
        q=ifftshift(q);
        out=zeros(size(in));  % allocate the output array
        for i=1:nim
            fq=q.*fftn(in(:,:,:,i));
            out(:,:,:,i)=real(ifftn(fq));
        end;

