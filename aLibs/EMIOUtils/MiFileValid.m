function ok=MiFileValid(name)

ok=true;
    
try
    if numel(name)>4 || isa(name,'struct')
        mi=name;
        name=[mi.infoPath mi.baseFilename 'mi.txt'];
    end;
    mi=ReadMiFile(name);
catch
    ok=false;
end;
