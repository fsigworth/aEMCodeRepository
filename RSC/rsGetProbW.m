function pW=rsGetProbW(ri,ialpha,ibeta,r0,yClick,xp)
% Compute the probability of the projected particle center w given alpha
% and beta values.

% ri.n=65;
% ialpha=-1;
% ibeta=8;
% r0=50;
% yClick=50;
% xp.sigmaR=0;
% xp.sigmaB=0;
% ri.angleSteps=[10 10 12];
% 
beta=ri.angleSteps(2)*(ibeta-1);
alpha=ri.angleSteps(1)*ialpha;

n=ri.n;

varX=(r0*pi/180*xp.sigmaR)^2+(r0*pi/180*ri.angleSteps(1))^2/12;
varY=(r0*pi/180*xp.sigmaR*abs(cosd(beta)))^2 ...
     +(xp.sigmaB*sind(beta))^2 ...
     +(r0*pi/180*ri.angleSteps(2))^2/12;

x0=r0*sind(beta)*sind(alpha);
y0=r0*sind(beta)*cosd(alpha)-yClick;

n0=-ceil((n-1)/2);
n1=n0+n-1;
[x y]=ndgrid(n0-x0:n1-x0,n0-y0:n1-y0);

pW=2/(2*pi*sqrt(varX+varY))...
    *exp(-(x.^2/(2*varX)+y.^2/(2*varY)));
imacs(pW);

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% sdFraction=.5;  % max spacing of projected pixels, relative to sigmaC
% boundarySDs=4;  % number of SDs of boundary in evaluation domain
% ds=4;  % oversampling of the probability function
% pf=2;  % padding factor
% 
% % Pick an angular step da so that r0*da < sdFraction*xParams.sigmaC
% da0=xp.sigmaC*r0*sdFraction*180/pi;  % maximum allowable step in degrees
% sts=ri.angleSteps(1:2)/2;   % magnitude of half-steps in alpha and beta
% da=sts./ceil(sts/da0);  % now an integer fraction of half-step
% 
% 
% 
% 
% 
% angleStep=ceil(minAngleStep*pi/180*r0/(xp.sigmaC));
% da0=xp.sigmaC/(r0*pi/90*minAngleStep);  % in fraction of steps
% if da0>1/3
%     da0=1/3;
% end;
% da=1/ceil(1/da0);  % integer fraction of 
% 
% else
%     da=floor(da0);
% end;
