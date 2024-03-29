function c=FRCorr(m1,m2,stackMode)
% c=FRCorr(m1,m2)
% Compute the fourier ring correlation between the images m1 and m2.
% The FRC is defined as
%           \sum{F1 .* conj(F2)}
% c(i) = -----------------------------
%        sqrt(\sum|F1|^2 * \sum|F2|^2)
%
% Where F1 and F2 are the Fourier components at the given spatial frequency i.
% The returned vector c has dimension n/2-1.
%
if nargin<3
    stackMode=0;
end;
n=size(m1);
nImgs=size(m2,3);

f1=fftshift(fft2(m1));
f2=fftshift(fft2(m2));

num=Radial(real(f2(:)'*f1(:)));
den=Radial(sqrt(f1(:)'*f1(:)+f2(:)'*f2(:)));

c=num./den;



d0=disco(n,0.5);
c=[];
for i=1:n/2-1
	d1=disco(n,i+0.5);
	ring=d1-d0;
	nr=sum(ring);
	r1=ring .* f1;
	r2=ring .* f2;
	num=real(sum(sum(r1.*conj(r2))));
	den=sqrt(sum(sum(abs(r1).^2))*sum(sum(abs(r2).^2)));
	c(i)=num/den;
	d0=d1;
end;
