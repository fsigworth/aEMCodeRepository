% MiWriteAll
% write out the allMis array to a set of mi files

outDir='Info_8/';
CheckAndMakeDir(outDir,1);

nmi=numel(allMis);
for i=1:nmi
    mi=allMis{i};
    miName=[outDir mi.baseFilename 'mi.txt'];
    disp([num2str(i) ' ' miName]);
    WriteMiFile(mi,miName);
end;
