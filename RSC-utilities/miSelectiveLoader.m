% miSelectiveLoader

checkMerged=1;

sourceDir='Info/';
targetDir='Info_m/';
suffix='m.mrc';

doExecute=1;
opString='cp -n';

disp(pwd)
% Load all the mi files from the target
miNames=f2FindInfoFiles(sourceDir);
miNames=sort(miNames);
nm=numel(miNames);
disp('Loading mi files...');
for i=1:nm
    %     disp(names{i});
    mis{i,1}=ReadMiFile(miNames{i});
end;
disp(' done.');
%%
CheckAndMakeDir([targetDir],1);
k=0;
allMis=cell(0);
for i=1:nm
    mi=mis{i};
    miName=[mi.baseFilename 'mi.txt'];
    mergedName=[mi.procPath mi.baseFilename suffix];
    if exist(mergedName,'file')
        k=k+1;
        allMis(k)=mis(i);
        str=[opString ' ' sourceDir miName ' ' targetDir miName];
        disp(str);
        if doExecute
            system(str);
        end;
    end;
end;
save([targetDir 'allMis.mat'],'allMis');
filesCopied=k
