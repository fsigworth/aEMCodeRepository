function [lastName,lastIter]=reFindLatestMoi(path,ri,iTwin,maxIter)
% To search for restarting data, we start with iteration 1 and search for
% an moi file with the highest iteration number, as constructed from
% reGetNameCode, assumed
% to be in the current directory.  If none is found, a null string is
% returned, along with iter=0.
% ---we should include a check for the other twin as well---

if nargin<4
    maxIter=150;
end;

path=AddSlash(path);
d=dir(path);
nd=numel(d);
names=cell(nd,1);
for i=1:nd
    names{i}=d(i).name;
end;
lastName='';
lastIter=0;
for iter=1:maxIter
    moiName=[reGetNameCode(ri,iter,iTwin) 'moi.mat'];
    if any(strcmp(moiName,names))
        lastIter=iter;
        lastName=moiName;
    end;
end;
