function [iAlpha iBeta]=rsGetAngleIndices(angles,nAlpha,nBeta)
% Create an array t x 2 of angle values which are ~uniform on a hemisphere,
% 
% If nBeta is odd, we include a point at 90 degees; if even, not.
% e.g. nBeta = 2 yields 0 and 60 degrees.
corr=(mod(nBeta,2)+1)/2; % .5 even, 1 odd
betaStep=90/(nBeta-corr); % odd nBeta gives a final value at 90 degrees;
% even nBeta gives a final value of 90*(nBeta-1)/(nBeta-.5))
alpha(alpha<0)=alpha(alpha<0)+360;
alpha(alpha>360)=alpha(alpha>360)-360;

iBeta=min(round(angles(2,:)/betaStep)+1,nBeta);
beta=betaStep*(iBeta-1);
na=min(nAlpha*round(sind(beta)*nAlpha+.5),nAlpha);
q=round(angles(1,:).*na/360);
q(q>na-.5)=0;  % assign the angles near 360 to 0.
iAlpha=min(round(q)+1,na);
