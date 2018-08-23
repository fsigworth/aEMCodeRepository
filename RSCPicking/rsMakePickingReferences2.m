% function [allTemplates membraneOffset angleInds ds]=rsMakePickingReferences2
% inds is the indices for the upper hemisphere alpha and beta only.

load /Volumes/TetraData/Structures/AMPAR/3KG2map5.8AGoodScale.mat
nt=size(map,1);
map=map*5.8;  % approx amplitude correction (V-A scaling)
membraneOffset=-24/2;  % downsampled map by 2.
ds=2;

% addNonParticles=0;
%%
nAlpha=32; % about 10 degrees
nBeta=12;  % even is best.  Number on hemisphere.
nGamma=8;
symmetry=2;
gammaStep=360/(symmetry*nGamma);

[hemiAngles angleInds]=rsListHemisphereAngles(nAlpha, nBeta);
nHemiAngles=size(hemiAngles,1);
% make two copies, one for upper and one for lower hemisphere.
sphereAngles=[hemiAngles
       repmat([0 180],nHemiAngles,1)-hemiAngles];
angleList=repmat(sphereAngles,nGamma,1);
    % increment gamma each 2*nHemiAngles
gammas=repmat((0:nGamma-1)*gammaStep,2*nHemiAngles,1);
angleList(:,3)=gammas(:);  % gamma angles are most slowly varying.
% angle list is of size
% (nHemiAngles x nH x nGamma, 3)  where nH=2 is the number of hemispheres.
% Note that the nHemiAngles goes from beta=0 to 89, then beta =  180 to 91;
% that is, the same projected position is described twice.  Then all this
% is repeated nGamma times, with the gamma angle incremented each time.

nAngles=size(angleList,1);
%%
disp(['Making ' num2str(nAngles) ' templates']);

tic
% allTemplates=rsMakeTemplatesQuick(angleList,map);
allTemplates=rsMakeTemplates(angleList,map);
toc
%%
angleList=reshape(angleList,nHemiAngles,2,nGamma,3);
allTemplates=reshape(allTemplates,nt,nt,nHemiAngles,2,nGamma);
% size(allTemplates)
% %%
% if addNonParticles
%     % put in the non-particles
%     nrds=size(npRefs,4);
%     allTemplates(:,:,:,nGamma+1:nGamma+nrds,1)=npRefs*100;
%     size(allTemplates)
%     % templates are indexed as (x,y,ptrs,gamma,mirror)
% else
%     nrds=0;
% end;
