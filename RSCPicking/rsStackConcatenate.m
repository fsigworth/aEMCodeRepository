function [accumSi, img2]=rsStackConcatenate(si, imgs, accumSi, img2)
% Concatenate the stackInfo structure si and the stack imgs into the
% accumulated structure and stack.  Do error checking to make sure the image
% sizes and numbers are consistent.
% 
[n,ny,nim]=size(imgs);
nmi=size(si.ctfs,3);

% Make sure the activeFlags are set up.
if ~isfield(si,'activeFlags')
    si.activeFlags=true(nim,1);
end;
if ~isfield(si,'activeFlagLog')
    si.activeFlagLog={date};
end;
if ~isfield(si,'miGroup') % labels each micrograph with a micrograph group number.
    si.miGroup=ones(nmi,1,'uint16');
end;

% if the accumulators are empty, just copy and return.
if numel(img2)==0
    img2=imgs;
    accumSi=si;
    return
end;

[na, nay, naccim]=size(img2);
nami=size(accumSi.ctfs,3);

if ~all([n ny]==[na nay])
    error('Inconsistent image sizes');
end;
if ~isfield(accumSi,'activeFlags')
    accumSi.activeFlags=true(naccim,1);
end;
if ~isfield(accumSi,'activeFlagLog')
    accumSi.activeFlagLog={date};
end;
if ~isfield(accumSi,'miGroup')
    accumSi.miGroup=ones(nami,1,'uint16');
end;
maxAccumMiGroup=max(accumSi.miGroup);

if abs(si.pixA-accumSi.pixA)>.01  % differences <.01 A are ok
    error(['Inconstent pixel sizes ' num2str([si.pixA accumSi.pixA])]);
end;

% Check for unique mis by looking at the identifier numbers
nmi=numel(si.mi);
nami=numel(accumSi.mi);
% idsa=zeros(nmi,1);
% for i=1:nmi
%     idsa(i)=accumSi.mi{i}.identifier;
% end;
% for i=1:numel(si.mi)
%     id=si.mi{i}.identifier;
%     if any(idsa==id)
%         error('Redundant mi structure found');
%     end;
% end;

% Test baseFilenames
aBaseNames=cell(nami,1);
for i=1:nami
    if numel(accumSi.mi{i})>0
        aBaseNames{i}=accumSi.mi{i}.baseFilename;
    else
        aBaseNames{i}='';
    end;
end;
for i=1:nmi
    if numel(si.mi{i})>0
    nm=si.mi{i}.baseFilename;
    if any(strcmp(nm,aBaseNames))
         error(['Redundant mi structure found for: ' nm]);
    end;
    end;
end;

% Copy in the mi list and the ctfs
accumSi.mi(nami+1:nami+nmi)=si.mi;
accumSi.ctfs(:,:,nami+1:nami+nmi)=si.ctfs;
accumSi.miGroup(nami+1:nami+nmi)=si.miGroup+maxAccumMiGroup;

% update the indices to point to the new locations
accumSi.miIndex(naccim+1:naccim+nim)=si.miIndex+nami;

% Copy all the fields that have one entry per image
accumSi.miParticle(naccim+1:naccim+nim)=si.miParticle;
accumSi.alpha0(naccim+1:naccim+nim)=si.alpha0;
accumSi.yClick(naccim+1:naccim+nim)=si.yClick;
accumSi.rVesicle(naccim+1:naccim+nim)=si.rVesicle;
accumSi.sVesicle(naccim+1:naccim+nim)=si.sVesicle;

% Concatenate the activeFlags.
% Make the activeFlags consistent in size by cloning the last column.
nat=size(si.activeFlags,2);
nact=size(accumSi.activeFlags,2);
if nat>nact
    accumSi.activeFlags(:,nact+1:nat)=repmat(accumSi.activeFlags(:,nact),1,nat-nact);
    accumSi.activeFlagLog=si.activeFlagLog;  % copy the longest log
elseif nact>nat
    si.activeFlags(:,nat+1:nact)=repmat(si.activeFlags(:,nat),1,nact-nat);
end;
accumSi.activeFlags(naccim+1:naccim+nim,:)=si.activeFlags;
% We remove any particle selection
accumSi.pastParticle=uint32(1:naccim+nim);
accumSi.origParticle=uint32(1:naccim+nim);
% Copy the images
img2(:,:,naccim+1:naccim+nim)=imgs;
