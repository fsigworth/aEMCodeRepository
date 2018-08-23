% EllipticalVesicleModel




% t=10;
% n=256;
% b=100;
% e=[-.3 0 0 0];
% 
% a=sqrt(1./(1-e.^2))*b;
% em=sqrt(1-(b-t).^2./(a-t).^2).*sign(e);
% ep=sqrt(1-(b+t).^2./(a+t).^2).*sign(e);
% 
% de=EllipsoidDensity(n,b+t,ep) - EllipsoidDensity(n,b-t,em);
% 
% imacs(de)




mi.pixA=1.7;

% membrane model
vLipid=1.6;
thk=60;

% thk=30;
rise=6;
rise=2;
% Create the model, which is sampled in units of the original pixel size.
nm0=ceil(30/mi.pixA)*2+1;  % array for vesicle model; 60A nominal
mi.vesicleModel=(fuzzymask(nm0,1,thk/mi.pixA/2,rise/mi.pixA)...
                -fuzzymask(nm0,1,thk/mi.pixA*0.3,rise/mi.pixA))...
    *vLipid*mi.pixA;  % units of V.A per voxel
a=.2;
for phi=0:.1:2*pi
v=EllipticalVesicleFromModel(512,100,[a 0 a*cos(phi) a*sin(phi)],mi.vesicleModel);
imacs(v);
drawnow;
end;