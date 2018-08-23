% ImageOverview.m

% rexp='.*ala\.mrc'
rexp='.*\.mrc';
names={};
j=0;
d=dir;
for i=1:numel(d)
    q=regexp(d(i).name,rexp);
    if ~d(i).isdir && numel(q)>0
        j=j+1;
        names{j}=d(i).name;
    end;
end;
numFiles=j;
imgs={};
for j=1:numFiles
    [m,pixA]=ReadEMFile(names{j});
    imgs{j,1}=m;
    imgs{j,2}=pixA;
end;

    