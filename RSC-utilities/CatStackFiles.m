% CatStackFiles
% put up file selectors repeatedly to get stack files, then write out the
% concatenated file.
fullNames={};
fileNames={};
siPath='';
str=' ';
while ~isnumeric(siPath)
    % Put up a file selector for files *si.mat,
    [fnames, siPath]=uigetfile('*si.mat',['Select' str 'si files'],'multiselect','on');
    if ~isnumeric(siPath)  % user hasn't clicked Cancel
        if ~iscell(fnames)
            fnames={fnames};
        end;
        for i=1:numel(fnames)
            newName=[AddSlash(siPath) fnames{i}];
            disp(newName);
            fullNames{end+1}=newName;
            fileNames{end+1}=fnames{i};
        end;
        basePath=ParsePath(siPath);
        cd(basePath);
    end;
    str=' more ';
end;
%%
opts=struct;
[si,imgs,allNames,siPath,altImgs]=reLoadStackFiles(opts,fullNames,'');

totalImages=size(imgs,3)
totalSiEntries=numel(si.miIndex)
%%
[outName,outPath]=uiputfile('*si.mat','Output si file name');
if ~isnumeric(outPath)
    [pa,nm,ex]=fileparts(outName);
    baseName=nm(1:end-3);
    siName=[AddSlash(outPath) baseName 'asi.mat'];
    disp(['Writing ' siName]);
    save(siName,'si');
    stackName=[AddSlash(outPath) baseName 'astack.mrc'];
    disp(['Writing ' stackName]);
    WriteMRC(imgs,si.pixA,stackName);
    altStackName=[AddSlash(outPath) baseName 'austack.mrc'];
    disp(['Writing ' altStackName]);
    WriteMRC(altImgs,si.pixA,altStackName);
    disp('done.');
end;