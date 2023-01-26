function miEraseLogs(names)
% Erases the log field in each of the given mi files. names is a cell array
% of strings, each something like Info/TheFile.txt
%
nLogFields=zeros(numel(names))-1;
disp(['Erasing log fields in ' num2str(numel(names)) ' mi files.'])
for i=1:numel(names)
    mi=ReadMiFile(names{i});
    nLogs=numel(mi.log);
    nLogFields(i)=nLogs;
    mi.log={};
    WriteMiFile(mi,names{i});
    if mod(i,100)==0
        disp(i);
    end;
end;
