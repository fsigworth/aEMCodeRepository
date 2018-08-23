function [dis,si,stack,markInds]=rspLoadActiveStacks(dis)
% function [dis,stack,markInds]=rspLoadActiveStacks(dis)
% For SimpleStackPruner, load stacks and (1) select the set to display, and
% (2) select the sets to label with marks.  Returned are the stack
% (interleaved if dis.showDualStacks is true) and the cell array of indices
% of particles in the stack to mark.

if dis.showDualStacks
    [si0, stack0, names, siPath, ustack0]=reLoadStackFiles;
    dis.stackStep=2;
else
    ustack0=[];
    [si0, stack0, names, siPath]=reLoadStackFiles;
    dis.stackStep=1;
end;
disp('loaded:');
disp(names);
j=size(si0.activeFlags,2);
if j>1 % possibly trim the activeFlags
    afIndex=PickActiveFlagSet(si0,'Display all particles from which selection? ');
    active=si0.activeFlags(:,afIndex);
    if dis.showDualStacks
        [si,stack,ustack]=rsStackSplit(si0.activeFlags(:,afIndex),si0,stack0,ustack0);
    else
        ustack=[];
        [si,stack]=rsStackSplit(si0.activeFlags(:,afIndex),si0,stack0);
    end;
else
    active=si0.activeFlags(:,1);
    si=si0;
    stack=stack0;
    ustack=ustack0;
end;
origIndices=find(active);  % indices of original stack images

afIndex=1;
k=0;
markInds=cell(0,1);
while afIndex>0
    afIndex=PickActiveFlagSet(si,'Pick a set of particles to mark (0 to quit)? ');
    if afIndex>0
        k=k+1;
        markInds{k,1}=origIndices(si.activeFlags(:,afIndex));
    end;
end;
nim=size(stack,3);
if dis.showDualStacks % we'll interleave them
    stack0=stack;
    usd=std(ustack(:)); % make the unsubtracted set dim
    stack(:,:,1:2:2*nim)=dis.umul*(ustack-(1-dis.umul)*2*usd);
    stack(:,:,2:2:2*nim)=stack0;
    dis.stackStride=2;
else
    dis.stackStride=1;
end;
end

function j=PickActiveFlagSet(si0,txt)
disp([num2str(size(si0.miIndex,1)) ' particles total.']);
j=size(si0.activeFlags,2);
for i=1:j
    if numel(si0.activeFlagLog)<i
        si0.activeFlagLog{i,1}='';
    end;
    disp([num2str([i sum(si0.activeFlags(:,i))]) '  ' si0.activeFlagLog{i}]);
end;
j=MyInput(txt,j);
end
