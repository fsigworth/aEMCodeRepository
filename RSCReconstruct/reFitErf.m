function [fit,pars]=reFitErf(y)
n=numel(y);
x=(0:n-1)';
y=y(:);
a0=[.5 n/2 n/5];
pars=Simplex('init',a0);
for i=1:100
    fit=(.5-pars(1))-(.5+pars(1))* erf( (x-pars(2))/pars(3) );
    err=y-fit;
%     plot([y f]);
%     drawnow;
    sse=err(:)'*err(:);
    pars=Simplex(sse);
%     [sse a]
end;

