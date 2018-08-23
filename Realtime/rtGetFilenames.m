function names=rtGetFilenames(path,pattern,givePath)
path=AddSlash(path);
if numel(path)==0
    d=dir;
else
    d=dir(path);
end;
if givePath
    ppath=path;
else
    ppath='';
end;
names={};
for i=3:numel(d)
    if any(strfind(d(i).name,pattern))
        names{end+1}=[ppath d(i).name];
    end;
end;
