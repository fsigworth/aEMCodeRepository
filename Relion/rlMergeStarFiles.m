% rlMergeStarFiles

nData=zeros(0,1);
names=cell(0,1);
blockData=cell(0,1);
%%
fi=0;  % file index
done=false;
while ~done
    disp(['Getting the name of file ' num2str(fi+1)]);
    [nm,pa]=uigetfile('*.star');
    if isnumeric(pa)
        break;
    end;
    cd(pa);
    fi=fi+1;
    % Find the size and fieldnames of each file's data
    disp(['Reading ' pa nm]);
    [blockNames,blockData(fi),ok]=ReadStarFile(nm);
    %%
    s0=blockData{fi};
    names{fi}=fieldnames(s0);
    nData(fi)=numel(s0.(names{fi}{1}));
end;
nSets=fi;

%%
% Find the common fieldnames
nNames=numel(names{1});
oks=true(nNames,1);
for i=1:nNames
    fieldName=names{1}{i};
    for j=1:nSets
        oks(i)=oks(i) && any(strcmp(fieldName,names{j}));
    end;
end;

commonNames=names{1}(oks);
nComNames=numel(commonNames);
disp([num2str(nComNames) ' common names.']);

% Copy the first file into the big struct.
s=struct;
for i=1:nComNames
    s.(commonNames{i})=s0.(commonNames{i});
end;
for j=2:nSets
    for i=1:nComNames
        s.(commonNames{i})=[s.(commonNames{i}); blockData{j}.(commonNames{i})];
    end;
end;

disp('Getting an output file name');
[outName,outPath]=uiputfile('*.star');
if isnumeric(outPath)
    disp('Writing canceled');
    return
end;

cd(outPath);
%%
disp(['Writing ' outPath outName]);
WriteStarFileStruct(s,blockNames{1}(6:end),outName);
disp('done.');
