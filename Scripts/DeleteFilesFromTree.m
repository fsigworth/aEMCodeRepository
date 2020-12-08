function DeleteFilesFromTree(currPath,nDirs)
if nDirs==0
    return
end;
nDirs=0;
pattern='.MATLAB*';
doExec=1;
d=dir(currPath);
for i=3:numel(d)
    if d(i).isdir && numel(d(i).name)>2
        nDirs=nDirs+1;
        newPath=[currPath AddSlash(d(i).name)];
        DeleteFilesFromTree(newPath,nDirs)
        str=['rm -v ' newPath pattern];
        if doExec
            system(str);
        else
                disp(str);
        end;
    end;
end;
