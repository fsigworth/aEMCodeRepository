function [sis,imgs,altImgs]=reLoadStackGroup(ri,si,iTwin,iGroup,activeFlags)
% Load the appropriate parts of a stack for a node.  The fields
% ri.stackPath and ri.stackName point to the original si file. \
% Then based on ri.nTwins, ri.nNodes and ri.nSlices, create cell arrays
% with nSlices elements of si structurs, image stacks and altImage stacks.
% The si.activeFlags(:,end) is used to select elements of the stack, unless
% the optional activeFlags argument is given.

if nargin>4
    activeFlags=si.activeFlags(:,end);
end;
activeInds=find(activeFlags);
twinInds=activeInds(iTwin:ri.nTwins:end);

nTPerNode=ceil(totalActive/ri.nNodes);
startInd=(iNode-1)*nTPerNode+iTwin;
endInd=min(iNode*nTPerNode,totalActive);
nodeInds=activeInds(startInd:ri.nTwins:endInd);




load([ri.stackPath ri.siName]);  % load the si structure

% Construct the names of the stack and alt stack .mrc files.
p=regexp(ri.siName,'si.mat');
if numel(p)<1
    error(['No match found in name ' siName]);
end;
stackName=[siName(1:p(end)-1) 'stack.mrc'];
altStackName=[siName(1:p(end)-1) 'ustack.mrc'];

startImg=nodeInds(1);
endImg=nodeInds(end);
shortFlags=false(endImg-startImg+1,1);
shortFlags(nodeInds)=true;

nim=endImg-startImg+1;
shortStack=ReadMRC([siPath stackName],startImg,nim);
if ri.useAltImgs
    altShortStack=ReadMRC([siPath altStackName],startImg,nim);
end;
nActiveImgs=sum(shortFlags);
twinIndices=shortIndices(iTwin:2:end);
indsPerSlice=ceil(numel(twinIndices)/ri.nNodes);

sis=cell(ri.nSlices,1);
imgs=cell(ri.nSlices,1);
altImgs=cell(ri.nSlices,1);

for iSlice=1:ri.nSlices
    endInd=min(iSlice*indsPerSlice,nActiveImgs);
    sliceInds=shortIndices((iSlice-1)*indsPerSlice+1:endInd);
    imgs{iSlice}=shortStack(:,:,sliceInds-startImg+1);
    if ri.useAltImgs
        altImgs{iSlice}=altShortStack(:,:,sliceInds-startImg+1);
    end;
    sliceFlags=false(numel(si.miIndex),1);
    sliceFlags(sliceInds)=true;
    sis{iSlice}=rsStackSplit(sliceFlags,si);
end;
