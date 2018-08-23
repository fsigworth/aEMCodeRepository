% miSelectiveLoader

% sourceDir='Kv_1b/';
sourceDir='/Volumes/D254.2/Temp/';
targetDir='Kv_1sel/';
infoPath=[targetDir 'Info/'];
mergedDir='Merged/';
targets={mergedDir mergedDir mergedDir};
% targets={mergedDir};
suffixes={'ms.mrc' 'mvs.mrc' 'rscc.mat'};

opString='cp -n';
doExecute=1;

% Load all the mi files from the target
miNames=f2FindInfoFiles(infoPath);
nm=numel(miNames);
disp('Loading mi files...');
for i=1:nm
%     disp(names{i});
    mis{i,1}=ReadMiFile(miNames{i});
end;
disp(' done.');
%%
for j=2:numel(targets)
        CheckAndMakeDir([targetDir targets{j}]);
    for i=1:nm
        mi=mis{i};
    name=[mi.baseFilename suffixes{j}];
    str=[opString ' ' sourceDir targets{j} name ' ' targetDir targets{j} name];
    disp(str);
    if doExecute
        system(str);
    end;
    end;
end;
    