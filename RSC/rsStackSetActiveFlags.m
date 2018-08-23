function si1=rsStackSetActiveFlags(activeFlags,si0,logString)
% function si1=rsStackSetActiveFlags(activeFlags,si0)
% If we have a boolean vector activeFlags derived by StackSplit of si0,
% which therefore only has elements for the si0 active images, we now
% expand the boolean to the proper size, add it as
% another column to si0.activeFlags and return it as si1.
if nargin<3
    logString='';
end;
si1=si0;
ai=size(si0.activeFlags,2);
oldActiveInds=find(si0.activeFlags(:,ai));
newActiveInds=oldActiveInds(activeFlags);
q=false(size(si0.miIndex));
q(newActiveInds)=true;
si1.activeFlags(:,ai+1)=q;
si1.activeFlagLog{ai+1,1}=[date '  ' logString];
end