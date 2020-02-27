% RenameFiles.m


inPattern='mv.mrc';
outPattern='mvs.mrc';

d=dir;
for i=1:numel(d)
    ptr=strfind(d(i).name,inPattern);
    if numel(ptr)>0
        outName=[d(i).name(1:ptr(end)-1) outPattern];
        disp(outName);
    end;
end;