function mi=ReadMiExtension(mi, miName)

if isfield(mi,'mieFields') && numel(mi.mieFields)>0 % An mie file exists
    [pa,nm]=fileparts(miName);
    name=[AddSlash(pa) nm 'e.mat'];
    if ~exist(name,'file')
        warning(['mie file not found: ' name]);
        return
    end;
    try
        load(name);
    catch
        warning(['mie file corrupted: ' name]);
        return
    end;
    if ~strcmp(mi.timestamp,mie.timestamp)
        warning(['Inconsistent mie timestamp: ' mie.timestamp{1}]);
    end;
    for i=1:numel(mi.mieFields)
        f=mi.mieFields{i};
        if isfield(mie,f)
            mi.(f)=mie.(f);
        end;
    end;
end;
