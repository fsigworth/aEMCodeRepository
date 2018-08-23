% function [me,sd,fit]=EMGaussFit(y,xvals,lims,ymax)
% Using the EM algorithm, Fit a Gaussian to part of a histogram.  y are the
% bin heights, xs are the
% bin centers, and lims gives the range over which y values are valid, i.e.
% we use y(lims(1):lims(2)) as the range to be fitted.  Returned are the
% mean and 

fc=.1;
q=GaussFilt(m0,fc);
[y,xs]=hist(q(:),400);
[ymax,p1]=max(y);
lims(1)=find(y>0.9*ymax,1,'last');
lims(2)=find(y>0.2*ymax,1,'last');
xvals=xs;
nargin=4;
% ymax=y(lims(1));
%
% n=256;
% p1=128;
% p2=256;
% xs=(1:256)';
% x1=1.0*exp(-((xs-60).^2/(2*25^2)));
% x2=1.0*exp(-((xs-100).^2)/(2*20^2));
% x=x1+x2;
%
% y0=zeros(n,1);
% y0(p1:p2)=x(p1:p2);
y=y(:);
xvals=xvals(:);
p1=lims(1);
p2=lims(2);

z0=zeros(size(y));
z0(p1:p2)=y(p1:p2);
me=xvals'*z0/sum(z0);
si=(z0'*(xvals-me).^2)/sum(z0);
z=z0;
f2=.5;
for i=1:1000
    newdist=exp(-(xvals-me).^2/(2*si));
    sumNew=sum(newdist(p1:p2));
    sumOld=sum(z(p1:p2));
    if nargin<4
        ymax=sumOld/sumNew;
    end;
    z=newdist*ymax;
    if i>f2
        f2=2*f2;
        bar(xvals,z0,'r'); hold on;
        plot(xvals,[z y]); hold off;
        title([me sqrt(si)]);
        drawnow;
    end;
    z(p1:p2)=y(p1:p2);
    me=xvals'*z/sum(z);
    si=sum((xvals-me).^2.*z)/sum(z);
    % pause
end;
fit=newdist*sumOld/sumNew;
sd=sqrt(si);
