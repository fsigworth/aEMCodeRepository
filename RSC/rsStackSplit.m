function [sis,imgss,altImgss]=rsStackSplit(activeFlags,si,imgs,altImgs,startImg)
% function [sis,imgss]=rsStackSplit(activeFlags,si,imgs)
% Takes a StackInfo structure and a stack of images (or a stack of stacks
% of rotated images) and selects a subset according to the activeFlags.
% The si.activeFlags field is not used in pruning, and sis.activeFlags is
% set to true(nim,1).
% The mi's and ctfs are pruned, but not not ctfGroups.  If the imgs
% argument is missing, only si is operated on.
% The imgs stack can be smaller than the entire stack described by si.
% If so, it is assumed to be a contiguous part of the entire stack,
% starting with startImg (default is 1).

nim=numel(si.miIndex);
if nargin>2
    sz=size(imgs);
    inputNim=sz(end);
end;    
if nargin<5
    startImg=1;
end;

% Prune the set of mi's
nmi=numel(si.mi);
w=ones(nim,1);
snim=sum(activeFlags);
if snim<1
    error('No active images');
end;
h=WeightedHisto(si.miIndex(activeFlags),w,nmi);
activeMis=h>0;
newMiIndices=cumsum(uint16(activeMis));

if ~isfield(si,'origParticle')
    si.pastParticle=(1:uint32(nim))';
    si.origParticle=si.pastParticle;
end;
if ~isfield(si,'origMi')
    si.origMi=uint16(1:nmi);
end;


sis=struct;

% Prune all the particle-enumerated arrays
sis.origParticle=si.origParticle(activeFlags);
sis.pastParticle=find(activeFlags);
sis.miIndex=newMiIndices(si.miIndex(activeFlags));
sis.miParticle=si.miParticle(activeFlags);
sis.alpha0=si.alpha0(activeFlags);
sis.yClick=si.yClick(activeFlags);
sis.rVesicle=si.rVesicle(activeFlags);
sis.sVesicle=si.sVesicle(activeFlags);
sis.activeFlags=si.activeFlags(activeFlags,:);  % modify the active flags too!
if isfield(si,'ctfGroupIndex')
    sis.ctfGroupIndex=si.ctfGroupIndex(activeFlags);
end;

% Copy the mi's and ctfs
sis.mi=si.mi(activeMis);
sis.origMi=si.origMi(activeMis);
sis.ctfs=si.ctfs(:,:,activeMis);

% Copy the invariant arrays
sis.pixA=si.pixA;
sis.weights=si.weights;
if isfield(si,'mbnOffset')
    sis.mbnOffset=si.mbnOffset;
end;
if isfield(si,'ctfGroups')
    sis.ctfGroups=si.ctfGroups;
end;

% Mark that we've modified the active flags field
sis.activeFlagLog=si.activeFlagLog;
sis.activeFlagLog(1,1)={['rsStackSplit  ' date]};
nafs=sum(si.activeFlags,1);
for i=2:numel(nafs)
    if nafs(i)>sum(sis.activeFlags(:,i))
        sis.activeFlagLog{i,1}=[si.activeFlagLog{i,1} ' (split)'];
    end;
end;

% Prune the stack, if that's desired

imgss=[];
altImgss=[];
if nargin>2 && nargout>1
    
    imgActiveFlags=activeFlags(startImg:startImg+inputNim-1);
    if sum(imgActiveFlags)<snim
        error('Inconsistent flags');
    end;

    % Copy the images
    switch numel(sz)
        case 3
            imgss=imgs(:,:,imgActiveFlags);
        case 4
            imgss=imgs(:,:,:,imgActiveFlags);
        otherwise
            error(['Wrong dimension for image stack: ' num2str(dimsIm)]);
    end;
end;
if nargin>3 && nargout>2 && numel(altImgs)>1
    % Copy the alt images
    switch numel(sz)
        case 3
            altImgss=altImgs(:,:,imgActiveFlags);
        case 4
            altImgss=altImgs(:,:,:,imgActiveFlags);
    end;
end;
