% FindBadInfos.m

allNames=FindInfoFiles;
doList=0;
modConst=20;
modConst2=500;
nmi=numel(allNames);
disp([num2str(nmi) ' Info files total.']);
badNames=cell(0,1);
j=0;
for i=1:nmi
    nm=allNames{i};
    mi=ReadMiFile(nm,'noiseModelPars');
%    mi=ReadMiFile(nm,'frameSets');
%     mi=ReadMiFile(nm,'weights');
%    if numel(mi.frameSets)<4 || mi.frameSets(1,2)<20
%     mergeOk=numel(mi.imageFilenames)>0 && exist([mi.procPath mi.baseFilename 'mvz.tif'],'file');
pwOk=isfield(mi,'noiseModelPars') && numel(mi.noiseModelPars)>0;
%    if numel(mi.frameSets)<2 || mi.frameSets(1,2)>100  || ~mergeOk % no jump
if ~pwOk
    %     if mi.weights(2)<1
    j=j+1;
    badNames{j,1}=nm;
    if doList
        disp([num2str(j) ' bad: ' nm]);
    else
        if mod(i,modConst)==0
            fprintf('x');
        end;
    end;
else
    if doList
        disp(['                                 ' nm]);
    elseif mod(i,modConst)==0
        fprintf('o');
    end;
end;
if ~doList && mod(i,modConst2)==0
    fprintf('\n%d: ',i);
end;
end;
fprintf('\n');
disp('done');
