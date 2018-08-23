% SimAngleCoverage
% Simulate the angular coverage of a sphere
% for Microscopy review paper.
multiAxis=1;
nphi=160;
ntheta=40;
if multiAxis
    refAngs=SphereAngles3(nphi,ntheta,1)*pi/180; % nr x 3
else
    % Make a single-axis rot
    dphi=180/nphi;
    refAngs=(0:nphi-1)'*dphi*pi/180;
    refAngs(1,3)=0;
    refAngs(:,2)=pi/2;
end;

nr=size(refAngs,1)
uz=[0 0 1];
refVecs=zeros(3,nr);
for i=1:nr
    refVecs(:,i)=(uz*EulerMatrixStandard(refAngs(i,:)))';
end;

nv=10000;
dpMin=zeros(nv,1);
for i=1:nv
    v0=randn(1,3);
    v0=v0/sqrt(v0*v0');  % unit 3d vector
    dps=v0*refVecs;
    dpMin(i)=min(abs(dps));
end;
hist(dpMin,100);
medianErr=median(dpMin);
meanErr=mean(dpMin);
percent90=Percentile(dpMin,.9);
percent99=Percentile(dpMin,.99);
maxErr=max(dpMin);

disp([meanErr percent90 percent99 maxErr]);