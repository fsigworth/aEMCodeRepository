% function [allTemplates membraneOffset inds ds]=rsMakePickingReferences

load /Volumes/TetraData/Structures/AMPAR/3KG2RotMap5.8A.mat
nt=size(map,1);
map=map/64;  % approx amplitude correction (V-A scaling)
membraneOffset=-24/2;  % downsampled map by 2.
ds=2;

addNonParticles=0;
%%
nAlpha=32; % about 10 degrees
nBeta=24;
nGamma=8;
symmetry=2;

[angleList inds]=rsListSphereAngles(nAlpha, nBeta);
nAngles=size(angleList,1);

% Add the gamma angles

dGamma=360/symmetry/nGamma;
allTemplates=zeros(nt,nt,nAngles,nGamma);
angles=angleList;
angles(1,3)=0;  % add 3rd dimension
disp(['Making ' num2str(nGamma*nAngles) ' templates']);
tic
    for i=1:nGamma
        gamma=(i-1)*dGamma;
        tangles=angles;
        tangles(:,3)=gamma;
        q=rsMakeTemplatesQuick(tangles,map);
        allTemplates(:,:,:,i)=q;
    end;
toc
% size(allTemplates)
%%
if addNonParticles
    % put in the non-particles
    nrds=size(npRefs,4);
    allTemplates(:,:,:,nGamma+1:nGamma+nrds,1)=npRefs*100;
    size(allTemplates)
    % templates are indexed as (x,y,ptrs,gamma,mirror)
else
    nrds=0;
end;
