function [gSi,gImgs,gAltImgs,groupSize,twinInds]=reLoadStackGroup(ri,si,iTwin,iGroup,useAltImgs)
% function [gSi,gImgs,gAltImgs,groupSize,twinInds]=reLoadStackGroup(ri,si,iTwin,iGroup,useAltImgs)
% Load, downsample and crop the appropriate parts of a stack for the group
% of images corresponding to iTwin and iGroup, i.e. all the images to be
% operated on by a single node. The si structure for the group is
% returned.  The fields ri.stackPath and ri.stackName are used to load the
% original si file. If altImgs are not available or requested, altImgs is 
% set to 0. The ri.twinFlags(:,iTwin) is used to select the
% relevant elements of the original stack, unless the optional activeFlags
% argument is given. The default for useAltImgs is true if the gAltImgs
% argument is given.
% The returned gImgs and gAltImgs are of size ri.nCurrent

if nargin<5
    useAltImgs=nargout>3;
end;

% Split off odd or even, for iTwin=1 or 2.
% (get every image if ri.nTwins=1 and iTwin=1)
% activeInds=find(activeFlags);
% twinInds=activeInds(iTwin:ri.nTwins:end);
twinInds=find(ri.twinFlags(:,iTwin));
nTwinInds=numel(twinInds);

% Get the Group for this node
groupSize=ceil(nTwinInds/ri.nGroups);
startInd=(iGroup-1)*groupSize+1;
endInd=min(iGroup*groupSize,nTwinInds); % the last group might be smaller.
groupInds=twinInds(startInd:endInd);
% nGroup=endInd-startInd+1;  % number of images in this group

% si=LoadStruct([ri.stackPath ri.siName]);  % load the si structure

% Construct the names of the stack and alt stack .mrc files.
p=regexp(ri.siName,'si.mat');
if numel(p)<1
    error(['No match found in name ' siName]);
end;
stackName=[ri.siName(1:p(end)-1) 'stack.mrc'];
altStackName=[ri.siName(1:p(end)-1) 'ustack.mrc'];

% Get ready to load only a subset of the entire stacks
startImg=groupInds(1);
endImg=groupInds(end);
nShort=endImg-startImg+1;

% Read the stack files
disp(['Reading ' ri.stackPath stackName]);
shortStack=ReadMRC([ri.stackPath stackName],startImg,nShort);
if useAltImgs
    shortAltStack=ReadMRC([ri.stackPath altStackName],startImg,nShort);
else
    shortAltStack=[];
end;

% Convert the group indices to flags
groupFlags=false(size(ri.twinFlags,1));
groupFlags(groupInds)=true;

% Get the gSi structure and substacks just for the group.
[gSi,gImgs,gAltImgs]=rsStackSplit(groupFlags,si,shortStack,shortAltStack,startImg);

% Downsample and crop the stacks
disp('Downsamping');
gSiA=gSi;
[gSi, gImgs]=rsStackCrop(gSi,gImgs,ri.nCropU);
[gSi, gImgs]=rsStackDownsample(gSi,gImgs,ri.nCurrent);
if numel(gAltImgs)>0
    [gSiA, gAltImgs]=rsStackCrop(gSiA,gAltImgs,ri.nCropU);
    [gSiA, gAltImgs]=rsStackDownsample(gSiA,gAltImgs,ri.nCurrent);
end;
disp('done.');
% % Make the flags for the slices
% nSlice=ceil(nGroup/ri.nSlices);  % number of images in a slice
% sliceFlags=false(nGroup,ri.nSlices);
% for i=1:ri.nSlices
%     startInd=(i-1)*nslice+1;
%     endInd=min(i*nSlice,nGroup);
%     sliceFlags(startInd:endInd,i)=true;
% end;

% if nargout>2
%     % Create the actual slices
%     sis=cell(ri.nSlices,1);
%     imgs=cell(ri.nSlices,1);
%     altImgs=cell(ri.nSlices,1);
%     for i=1:ri.nSlices
%         [sis{i},imgs{i},altImgs{i}]=rsStackSplit(sliceFlags(:,i),groupSi,gImgsC,gAltImgsC);
%     end;
% end;