function [msk varMap]=AutoHoleMask(m,pixA,pars,varMode,varMap)
% pars.n0
% pars.thresh
% pars.width
% pars.edge

n=size(m);
n0=pars.n0;
binFactor=round(n(1)/n0);
frac=1.2*pars.thresh;
nw=n/binFactor;
width=n0/10*pars.width;
erw=width;
edge=n0/4*pars.edge;
fcb=.05;

if varMode
    if nargin<5  % varMap not given, so let's compute it.
        fm=fftn(m);
        hexp=8;
        fh=.005;
        fl=pixA/30;
        f=RadiusNorm(n)/pixA;  % Frequency in inverse A
        % Butterworth bandpass
        H=1./(1+(fh./f).^hexp+(f./fl).^hexp);
        % Compute the local variance, binned
        varMap=BinImage(abs(ifftn(fm.*ifftshift(H)).^2),binFactor);
    end;
else
    varMap=GaussFilt(BinImage(-m,binFactor),fcb);
end;

q=varMap(:);
nq=numel(q);
% h=hist(q.^.2,1000);
% semilogy(h);
qsrt=sort(q);
thrMin=qsrt(round(nq*.25));  % get the 1% value
thrMax=qsrt(round(nq*.92));  % get the 99% value
% thr=qsrt(round(nq*frac));
thr=thrMin+frac*(thrMax-thrMin);

msk=varMap>thr;
expmask=GaussFiltDCT(msk,.133/width)>.1;
em2=GaussFiltDCT(expmask,.133/erw)>.9;

%% Extend to edge

emsk=1-em2;   % =1 in interior
mr=emsk;
for irot=1:4
    for j=1:nw(2)  % scan each row
        if mr(1,j)==1  % a border pixel is enabled
            q=find(mr(:,j)==0,1);
            if numel(q)>0 && q<edge
                mr(1:q,j)=0;
            end;
        end;
    end;
    mr=rot90(mr);
end;
msk=mr;
