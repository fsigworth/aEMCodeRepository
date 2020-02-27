% miEPURenamer.m
% Given a set of mi files, create a new set with simplified names, e.g.
% GridSquare_28359688_FoilHole_28373575_01mi.txt is changed to
% g688_h75_e01_mi.txt

nmi=numel(allMis);
strs=repmat("",1,5);
for i=1:nmi
 strs(i,:)=strsplit(string(allMis{i}.baseFilename),'_');
end;
strs(:,[1 3])=[]; % get rid of text fields
% strs=string(A);
lens=max(strlength(strs));

% find the minimum unique lengths
startingChar=zeros(1,3);


for j=1:3;
    cArray=char(strs(:,j));
    disp(lens(j));
    for i=1:lens(j)
        disp([num2str([j i]) '  ' (unique(cArray(:,i)))']);
        if numel(unique(cArray(:,i)))>1
            startingChar(j)=i;
            break
        end;
    end;
    newCArray=cArray(:,i:end);
    
%     i
end;


return




% miSelectiveLoader

checkMerged=1;

sourceDir='Info_m/';
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
