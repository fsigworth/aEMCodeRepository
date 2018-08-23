% EllipticalVesicleModel

t=10;
n=256;
b=100;
e=.2;

a=sqrt(1/(1-e^2))*b;
em=sqrt(1-(b-t)^2/(a-t)^2);
ep=sqrt(1-(b+t)^2/(a+t)^2);

de=EllipsoidDensity(n,b+t,ep) - EllipsoidDensity(n,b-t,em);

imacs(de)
