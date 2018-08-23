function [angles,isos]=rsGetRandomAngles(na,angleLimits,pIso)
% Create random angles of size na x 3 for particle-image simulation
% e.g. by calling reMakeTemplates(vol,angles).
% All angles are in degrees.
% the matrix angles consists of column vectors [alpha beta gamma] where
% - alpha is uniformly distributed from -angleLimits(1) to angleLimits(1),
% with 180 added with probability pIso
% - beta is sine distributed from angleLimits(2) to 180-angleLimits(2)
% - gamma is uniformly distributed from 0 to angleLimits(3)
% The angles are sorted in increasing beta, first all ~iso and then iso.

isos0=(rand(na,1)<pIso);
alphas0=180*isos0+angleLimits(1)*(2*rand(na,1)-1);
betas0=180/pi*acos((2*rand(na,1)-1)*cosd(angleLimits(2)));  % sin distributed
[~,inds]=sort(betas0+1000*isos0);  % sort isos first, then betas
betas=betas0(inds);
alphas=alphas0(inds);
isos=isos0(inds);
gammas=angleLimits(3)*rand(na,1);
angles=single([alphas betas gammas]);
