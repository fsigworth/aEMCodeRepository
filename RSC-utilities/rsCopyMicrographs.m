% rsCopyMicrographs
% From a set of mi files, find micrograph names and copy them from one directory to
% another.
inDir='/Volumes/D254/180226/Kv_1/Micrograph/'
outDir='/Volumes/D254/180226/Kv_1sel/Micrograph/'
imIndices=[1 2];
doExec=1;
op='mv';
maxEntries=inf;
d=dir('Info/');
nPicks=0;
nEntries=min(numel(d),maxEntries);
startEntry=1;
nmi=0;
miPicks=zeros(nEntries,1);
miDef=zeros(nEntries,1);
for i=startEntry:nEntries    
   name=['Info/' d(i).name];
   if strndcmp(name,'mi.txt')
       nmi=nmi+1;
       mi=ReadMiFile(name);
       for j=1:numel(imIndices)
           index=imIndices(j);
           mcName=mi.imageFilenames{index};
           str=[op ' ' inDir mcName ' ' outDir mcName];
           disp(str);
           if doExec
               system(str);
           end;
       end;
  end;
end;
       