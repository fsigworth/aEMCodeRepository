function [accumSi, accumImgs]=rsStackConcatenate(si, imgs, accumSi, accumImgs)
% function [accumSi, accumImgs]=rsStackConcatenate(si, imgs, accumSi, accumImgs)
% Concatenate the stackInfo structure si and the stack imgs into the
% accumulated structure and stack.  Do error checking to make sure the image
% sizes and numbers are consistent.
% 
[n,ny,nim]=size(imgs);

% if the accumulators are empty, just copy and return.
if numel(accumImgs)==0
    accumImgs=imgs;
    accumSi=si;
    return
end;

[na, nay, naccim]=size(accumImgs);
if ~all([n ny]==[na nay])
    error('Inconsistent image sizes');
end;

if si.pixA ~= accumSi.pixA
    error(['Inconstent pixel sizes ' num2str([si.pixA accumSi.pixA])]);
end;

% Check for unique mis by looking at the identifier numbers
nmi=numel(accumSi.mi);
idsa=zeros(nmi,1);
for i=1:nmi
    idsa(i)=accumSi.mi{i}.identifier;
end;
for i=1:numel(si.mi)
    id=si.mi{i}.identifier;
    if any(idsa==id)
        error('Redundant mi structure found');
    end;
end;

% Copy in the mi list and the ctfs
accumSi.mi(nami+1:nami+nmi)=si.mi;
accumSi.ctfs(:,:,nami+1:nami+nmi)=si.ctfs;

% update the indices to point to the new locations
si.miIndex=si.miIndex+nami;

% Copy all the fields that have one entry per image
accumSi.miIndex(naccim+1:naccim+nim)=si.miIndex;
accumSi.miParticle(:,:,naccim+1:naccim+nim)=si.miParticle;
accumSi.alpha0(:,:,naccim+1:naccim+nim)=si.alpha0;
accumSi.yClick(:,:,naccim+1:naccim+nim)=si.yClick;
accumSi.rVesicle(:,:,naccim+1:naccim+nim)=si.rVesicle;
accumSi.sVesicle(:,:,naccim+1:naccim+nim)=si.sVesicle;

accumSi.origMi=(1:numel(accumSi.mi))';
accumSi.origParticle=(1:numel(accumsi.miIndex))';
accumSi.pastParticle=accumSi.origParticle;

% Copy the images
accumImgs(:,:,naccim+1:naccim+nim)=imgs;
