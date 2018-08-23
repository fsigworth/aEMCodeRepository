function flags=reGetSliceFlags(ri,nGroup,iSlice)
% Get the binary flags for slices from a group of nGroup images. The
% returned flags is a logical vector nGroup x 1 in size. if iSlice==0, the
% flags for all the groups are returned in a logical array 
% nGroup x ri.nGroups in size.

% Handle the case for returning all the flags.
if nargin<3 || iSlice==0
    ng=ri.nSlices;
    is=1;
else  % return flags for one slice
    ng=1;
    is=iSlice;
end;
flags=false(nGroup,ng);
sliceSize=ceil(nGroup/ri.nSlices);
for ind=1:ng
    ig=ind-1+is;
    startInd=(ig-1)*sliceSize+1;
    endInd=min(ig*sliceSize,nGroup); % the last group might be smaller.
    sliceInds=startInd:endInd;
    flags(sliceInds,ind)=true;
end;
