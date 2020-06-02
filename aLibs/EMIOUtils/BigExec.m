function BigExec(cmd,sourceDir,endPattern,targetDir)
% Execute a linux command for a subset of files. Example:
% BigExec('mv','Merged','ms.mrc','Merged_s') is equivalent to
% mv Merged/*ms.mrc Merged_s/
% but has no limit on argument sizes etc.

d=dir(sourceDir);
d=d(1:20);
for i=1:numel(d)
    p=d(i);
    if ~d(i).isdir
        nm=d(i).name;
        if strndcmp(nm,endPattern,numel(endPattern))
            str=[cmd ' ' AddSlash(sourceDir) nm ' ' AddSlash(targetDir)];
            system(str);
        end;
    end;
end;
