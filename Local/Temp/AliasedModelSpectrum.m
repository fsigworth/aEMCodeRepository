function sp=AliasedModelSpectrum(n,f0,ex)
na=1;
ct=ceil(n+1)/2;
sp=zeros(n,n);
for i=1:2*na+1
    x=ct+(i-1-na)*n;
    for j=1:2*na+1
        y=ct+(j-1-na)*n;
        f=RadiusNorm(n,[x,y]);
        sp=sp+exp(-f/f0)+ex;
    end;
end;

