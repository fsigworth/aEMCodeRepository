% MiLoadAll.m
% Use the file selector to choose a set of mi files.  Load all of them into
% a single cell array, and store all their names too.
getNamesOnly=0;

[names, pathName]=uigetfile({'*mi.txt';'*mi.mat'},'Select mi files','multiselect','on');
if isnumeric(pathName) % File selection was cancelled
    return
end;
if ~iscell(names)
    names={names};
end;
[rootPath, infoPath]=ParsePath(pathName);
cd(rootPath);
%%
allMis=cell(0);
nNames=numel(names);
disp([num2str(nNames) ' mi files.']);
allNames=cell(nNames,1);
for i=1:nNames
    disp(names{i});
%     allNames{i}=[infoPath names{i}];
    allNames{i}=[infoPath names{i}];
    if ~getNamesOnly
        allMis{i,1}=ReadMiFile(allNames{i});
    end;
end;
save([infoPath 'allNames.mat'],'allNames');
disp('allNames.mat saved.');
if ~getNamesOnly
    save([infoPath 'allMis.mat'],'allMis','allNames');
    disp('allMis.mat saved.');
end;

