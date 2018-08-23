% ConvertInfoFiles

cd /Users/fred/EMWork/Hideki/170609
inDir='Info/';
outDir='Info_ext/';

di=dir(inDir);
CheckAndMakeDir(outDir,1);
tic
for i=199:numel(di)
    name=di(i).name;
    if any(strfind(name,'mi.txt'))
        disp(name);
        mi=ReadMiFile([inDir name]);
        WriteMiFile(mi,[outDir name]);
    end;
end;
toc
        