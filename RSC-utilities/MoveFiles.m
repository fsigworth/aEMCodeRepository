% MoveFiles.m
doExec=1;
nChars=12;
suffix='*m.mrc';

operator='cp';
% cd /ysm-gpfs/scratch60/fjs2/170417/KvLipo134_4/
cd ~/project/180226/

sourceDir='Kv_1/Merged/';
targetDir='Kv_1sel/Merged/';
targetBaseDir='Kv_1sel/';

% pick up the info file names in the target
load([targetBaseDir 'allNames.mat']);
% get rid of the 'Info/' and remove extension.
j=0;
for i=1:numel(allNames)
    [pa,nm,ex]=fileparts(allNames{i});
    if strcmp(ex,'.txt')
        j=j+1;
        infoNames{j}=nm;
    end;
end;
nTargetNames=j

%sourceDir='sq02w11/Micrograph_bad/';
ds=dir(sourceDir);
% %% Copy according to source dir entries.
% for i=1:numel(ds)
%     % check for a match of the first nc characters of the name
%     sourceName=ds(i).name;
%     q=strncmp(sourceName,targetNames,nChars); % see if a source file matches a target name.
%     if any(q)
%         name=[sourceName(1:nChars) suffix];
%         str=['!' operator ' ' sourceDir name ' ' targetDir];
%         disp(str);
%         if doExec
%             eval(str);
%         end;
%     end;
% end;
% return
% 

%% Copy according to target, no source
sourceNames=cell(numel(ds),1);
for i=1:numel(ds)
    sourceNames{i}=ds(i).name;
end;

for i=1:nTargetNames
    % check for a match of the first nc characters of the name
    
    q=strncmp(infoNames{i},sourceNames,nChars); % see if a source file matches a target name.
    name=[infoNames{i}(1:nChars) suffix];
    if any(q)
        str=['!' operator ' ' sourceDir name ' ' targetDir];
        disp(str);
        if doExec
            eval(str);
        end;
    else
        disp(['Not found: ' name]);
    end;
end;

return









for i=1:numel(nm)
str=['!mv ' sourceDir nm{i} ' ' targetDir nm{i}];
disp(str);
if doExec
    eval(str);
end;

end;
