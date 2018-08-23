% MoveFiles.m
doExec=1;
nChars=12;

cd /ysm-gpfs/scratch60/fjs2/170417/KvLipo134_4/
load sq02Names.mat
nm=sq02Names;
nNames=numel(nm);

sourceDir='sq02w11/Micrograph_bad/';
ds=dir(sourceDir);
targetDir='sq02w11/Micrograph/';

for i=1:numel(ds)
    % check for a match of the first nc characters of the name
    fName=ds(i).name;
    q=strncmp(fName,nm,nChars);
    if any(q)
        str=['!mv ' sourceDir fName ' ' targetDir fName];
        disp(str);
        if doExec
            eval(str);
        end;
    end;
end;
return
        










for i=1:numel(nm)
str=['!mv ' sourceDir nm{i} ' ' targetDir nm{i}];
disp(str);
if doExec
    eval(str);
end;

end;
