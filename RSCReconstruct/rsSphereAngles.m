function angles=rsSphereAngles(nAngles,symmetry,minBeta,gammaStepFactor)
% function angles=rsSphereAngles(nAngles,symmetry,minBeta,gammaStepFactor)
% create angles tiling a band about the equator of a sphere.
% the returned angels are of the form nAngles x [0 beta gamma].
% Betas are sampled more sparsely away from the equator, while gamma values
% are uniformly spaced.

if nargin<2
     minBeta=30;
end;
if nargin<3
     gammaStepFactor=1;
end;

maxGamma=360/symmetry;
% beta values are spaced as 1/sin(beta)
ys=zeros(nAngles,1);
y=0;
for i=1:nAngles
    b=minBeta+(i-1)/nAngles*(180-2*minBeta);
    y=y+1/sind(b);   
    ys(i)=y;
end;
dBeta=(180-2*minBeta)/y; % rise in beta per gamma step at beta=90
betas=minBeta+dBeta*ys;  % actual betas, widely spaced ner min and max
totalArea=dBeta*nAngles*maxGamma/gammaStepFactor;
qBeta=sqrt(totalArea/nAngles);  % size of box in beta direction
dGamma=qBeta*gammaStepFactor;   % gamma direction: times aspect ratio

angles=zeros(nAngles,3);
angles(:,2)=minBeta+dBeta*ys;
angles(:,3)=mod((0:nAngles-1)'*dGamma,maxGamma);

% plot(angles(:,3),angles(:,2),'.-');