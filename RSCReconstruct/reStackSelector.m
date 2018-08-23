% StackSelector
% assumes these variableds
% fname  -- cell array having si name.
% cls
% si
% allImgs

nImgs=numel(si.miIndex);
flags=false(nImgs,1);
for i=1:numel(goodClasses)
    q=cls{i};
    for j=1:numel(q)
        flags(q(j))=true;
    end;
end;
si0=si;
[si,imgsEd]=rsStackSplit(flags,si0,allImgs);
siName=['Ed' fname{1}]
save(siName,'si');
si=si0;
        p=regexp(siName,'si.mat');
        stackName=[siName(1:p(end)-1) 'stack.mrc'];
        stackName
       WriteMRC(imgsEd,si.pixA,stackName);
       
 
    