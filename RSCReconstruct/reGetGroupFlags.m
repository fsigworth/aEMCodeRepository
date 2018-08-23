function flags=reGetGroupFlags(ri,iTwin,iGroup)
% Get the binary flags for group collection of images, a logical vector
% nTwin x 1 in size.  nTwin is the number of images in the already squeezed stack.
% if iGroup==0, the flags for all the groups are returned in a logical
% array nTwin x ri.nGroups in size.

% Handle the case for returning all the flags.
if nargin<3
    iGroup=0;
end;
if iGroup==0
    ng=ri.nGroups;
    iGroup=1;
else
    ng=1;
end;
nt=ri.nTwin(iTwin);
flags=false(nt,ng);
groupSize=ceil(nt/ri.nGroups);
for ind=1:ng
    ig=ind-1+iGroup;
    startInd=(ig-1)*groupSize+1;
    endInd=min(ig*groupSize,nt); % the last group might be smaller.
    groupInds=startInd:endInd;
    flags(groupInds,ind)=true;
end;
