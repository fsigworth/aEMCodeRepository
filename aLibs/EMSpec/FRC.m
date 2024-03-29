function c=FRC(m1,m2)
% function c=FRC(m1,m2)
% Compute the fourier ring correlation between the images m1 and m2.
% The FRC is defined as
%           \sum{F1 .* conj(F2)}
% c(i) = -----------------------------
%        sqrt(\sum|F1|^2 * \sum|F2|^2)
%
% Where F1 and F2 are the Fourier components at the given spatial frequency
% i, and sums are taken over rings of radius i-1 <= r < i
% The returned vector c has dimension floor(n/2).
% If m1 and m2 have 3 dimensions, they are taken to be stacks of images.
% Then each sum is also taken
% over the stack of images, so a power-weighted FRC is obtained.
%
sz1=size(m1);
sz2=size(m2);
if numel(sz1)~=numel(sz2) || any(sz1~=sz2)
    error(['Unmatched dimensions of m1, m2: ' num2str(sz1) ',  ' num2str(sz2)]);
end;

f1=fftshift(fft2(m1));
f2=fftshift(fft2(m2));

num=real(Radial(sum(f1.*conj(f2),3)));
d1=Radial(sum(f1.*conj(f1),3));
d2=Radial(sum(f2.*conj(f2),3));
c=num./(sqrt(d1.*d2));

%
%
% d0=disco(n,0.5);
% c=[];
% for i=1:n/2-1
% 	d1=disco(n,i+0.5);
% 	ring=d1-d0;
% 	nr=sum(ring);
% 	r1=ring .* f1;
% 	r2=ring .* f2;
% 	num=real(sum(sum(r1.*conj(r2))));
% 	den=sqrt(sum(sum(abs(r1).^2))*sum(sum(abs(r2).^2)));
% 	c(i)=num/den;
% 	d0=d1;
% end;
