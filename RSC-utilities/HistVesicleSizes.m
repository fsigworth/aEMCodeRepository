% ShowVesicleSizes

allNames=nmsY22;
rvals=0;
for i=1:numel(allNames)
    nm=allNames{i};
%     disp(nm);
    mi=ReadMiFile(nm);
    if isfield(mi,'vesicle') && isfield(mi.vesicle,'r') & numel(mi.vesicle.r)>0
        rvals=[rvals;mi.vesicle.r(:,1)];
    end;
end;
numel(rvals)
rvals(rvals<0)=NaN;
rvals(rvals>500)=NaN;
hist(rvals,200);
