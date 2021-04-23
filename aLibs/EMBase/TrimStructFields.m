function sOut=TrimStructFields(sIn,iStart,iEnd);

fields=fieldnames(sIn);
sOut=struct;
for i=1:numel(fields)
    name=fields{i};
    sOut.(name)=sIn.(name)(iStart:iEnd,1);
end;
