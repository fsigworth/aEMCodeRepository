function out=BFilter(in, B, fc)
m=size(in);
ndim=sum(m>1); % number of non-singleton dimensions

% The unit of frequency in the various dimensions will be
% 1/m(1), 1/m(2), etc.  We want the output to be 1/sqrt(2) at fc, i.e.
% at fx=fc*m units.  The output will be exp(-(x.^2/fx.^2)*ln(2)/2);

b=B/4;
k=-B/(4*fc^2*m.^2);

k=1.782./fw;
n=size(in,1);
[x,y,z]=ndgrid(-n/2:n/2-1)
r=fftshift(RadiusNorm(m));

h=0.5*exp(-b*r.^2).*(1-erf(k*(r-fc)));
            fq=h.*fftn(in);
            if isreal(in)
                out(:,:,:,i)=real(ifftn(fq));
            else
                out(:,:,:,i)=ifftn(fq);
            end;
        end;
    end;
    
end;
