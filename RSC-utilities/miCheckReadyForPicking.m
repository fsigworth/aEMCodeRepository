% rsCheckReadyForPicking.m

load allNames
badNames={};
for i=1:numel(allNames)
    mi=ReadMiFile(allNames{i});
    inds=miDecodeLog(mi);
    if inds(6)<inds(5) || inds(5)<1 % not prewhitened, or no vesicle refinement
        badNames(end+1)=allNames(i);
        disp(allNames{i});
    end;
end;
