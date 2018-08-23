% CentroidPrecision

nv=10;
offs=.01;
vals=zeros(nv,1);
np=1e9;
sigma=8;
for i=1:nv
z=mean(round(sigma*randn(np,1)+offs));
vals(i)=z;
disp(z);
end;

MeanAndSD=[mean(vals) std(vals)]
