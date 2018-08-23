function [newNames,allNames]=FindNewFiles(oldFNames,extensions)
if nargin<1
    oldFNames=[];
end
if nargin<2
    extensions={'./'};
end;

% oldFNames=origFNames;
% extensions={'.mrc'};

pars.displayOn=1;
pars.extensions=extensions;
% startDir='./';

allNames=SearchDirectories(startDir,{},pars);

nf=numel(allNames);
newNames={};
nNew=0;
for i=1:nf
    testName=allNames{i};
    if ~any(strcmp(testName,oldFNames))
        newNames{end+1,1}=testName;
    end;
end;

