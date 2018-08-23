function [R,freqs]=RadialCTF(w,P,pixA)
% function [R,freqs]=RadialCTF(w,P,pixA)
% Given the CTF parameters, compute the appropriate 1D representation of
% the n x n, zero-center-frequency ctf or power spectrum w.  When there is no astigmatism, 
% this gives the same result as Radial(w).
% R(1) corresponds to frequencies [0,.5) *df; R(2) is frequencies [.5,1.5)*df and so
% on.  freqs are the corresponding bin centers 0, 1, 2... * df, where
% df=n*pixA.
% We use eqn. (8) of Penczek et al, Ultramicroscopy 2014.

% % Test code
% P=P0;
% w=c.^2;
% pixA=2.5;
% P.deltadef=0;
% P.theta=pi*rand;
% c=CTF(1024,pixA,P);
% w=c.^2;
%
% subplot(231);
% imacs(w);
% title(theta);

n = size(w,1);
org=ceil((n+1)/2);
[X,Y]=ndgrid(1-org:n-org,1-org:n-org);  % Make zero at x,y
r2=(X.^2+Y.^2);
alpha=atan2(Y,X);
pixA=double(pixA);  % everything has to be double else we have numerical problems.
S2=r2/(n*pixA)^2;

alpha0=-pi/4-double(P.theta);
d=double(P.defocus);
deltaD=2*double(P.deltadef);
cs=double(P.Cs);
la=double(P.lambda);

u=n*pixA/(cs*la)*sqrt(cs*(d-sqrt(d^2+cs^2*la^4*S2.^2 ...
    -2*cs*la^2*S2.*(d-deltaD/2*sin(2*(alpha+alpha0))))));
if any(imag(u(:))~=0)
    warning('WtHisto: u is complex');
    u=abs(u);
end;

[R, norm]=WeightedHisto(floor(u+1.5),w,floor(n/2));
R=R./(norm+eps);
freqs=(0:floor(n/2)-1)'/(n*pixA);
% subplot(232);
% plot([R Radial(w)]);
% drawnow;
