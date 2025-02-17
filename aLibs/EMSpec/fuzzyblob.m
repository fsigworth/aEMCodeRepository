function m = fuzzyblob(n,c,risetime)
% function m = fuzzyblob(n,a,risetime,origin)
% Create a disc of radius c(1), made with an error
% function with effective risetime tr.  The disc is distorted according to
% the Fourier series with complex coefficients c(1)=radius, c(2)=shift from
% center, c(3) ellipticity etc.
% The old fuzzydisc has an fc argument that corresponds to tr in this way:
% fc=0.1 <--> tr=3.3
% fc-0.2 <--> tr=1.7
% The origin variable is a vector; its default value is n/2+1 in each
% element.
% if r0 is a vector, an elliptical disc is made.  r0(1) is the basic
% radius; r0(2) is the x/y ellipticity; r0(3) is the diagonal ellipticity.
% r0=[rr 1 0] is circular.

ctr=n/2+1;  % default center is appropriate for ffts.
k=1.782/risetime;
org=[ctr+real(c(2)) ctr+imag(c(2))];
%             origin = [real(a(2))+ctr imag(a(2))+ctr];  % Coordinate of the center
%          z=(x-1i*y)./max(1,sqrt(x.^2+y.^2));  % complex angle
[r t]=Radius(n,org);
z=(exp(-1i*t));
a=single(zeros(n,n));
for i=1:numel(c)
    if i~=2
    a=a+real(c(i)*z.^(i-1));
    end;
end;
m=0.5*(1-erf(k*(r-a)));
