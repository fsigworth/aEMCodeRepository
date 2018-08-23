function flags=reGetFlags(ri,nActive,iTwin,iGroup,iSlice)
% Get the binary flags for the twin, group or slice collection of images.
% nActive is the number of images in the already squeezed stack.
% For twin, if iGroup==0;  for group if iSlice==0, etc.

% Figure out the twin flags
flags=false(nActive,1);
flags(iTwin:ri.nTwins:end)=true;
if iGroup==0  % we don't want a set of group flags
    return
else
    % Get the Group for this node
    twinInds=find(flags);
    nTwinInds=numel(twinInds);
    groupSize=ceil(nTwinInds/ri.nGroups);
    startInd=(iGroup-1)*groupSize+1;
    endInd=min(iGroup*groupSize,nTwinInds); % the last group might be smaller.
    groupInds=twinInds(startInd:endInd);
    flags=false(nActive,1);
    flags(groupInds)=true;
    if iSlice==0
        return
    else
        nGroup=numel(groupInds);
nSlice=ceil(nGroup/ri.nSlices);  % number of images in a slice
flags=false(nActive,1);
=false(nGroup,ri.nSlices);
for i=1:ri.nSlices
    startInd=(i-1)*nslice+1;
    endInd=min(i*nSlice,nGroup);
    sliceFlags(startInd:endInd,i)=true;
        
    end;
% nGroup=endInd-startInd+1;  % number of images in this group

nGroup=size(gImgs,3);
nSlice=ceil(nGroup/ri.nSlices);  % number of images in a slice
sliceFlags=false(nGroup,ri.nSlices);
for i=1:ri.nSlices
    startInd=(i-1)*nslice+1;
    endInd=min(i*nSlice,nGroup);
    sliceFlags(startInd:endInd,i)=true;
end;
