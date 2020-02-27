% miSortByDefocus

%  We assume we have allNames and allMis loaded, e.g. by running MiLoadAll,
% and we're in the expt directory.
% We create directories like Info__0.5 Info__1.0 etc. and copy info files
% into them.
% % MiLoadAll;

%%
fracMin=0.1; % Merge bins having a smaller fraction than this of the mode
doWrite=1;

nmi=numel(allMis);
defs=zeros(nmi,1);
for i=1:nmi
    defs(i)=allMis{i}.ctf(1).defocus;
end;
maxDef=max(defs);
minDef=min(defs);


quant=.5;
binMaxs0=ceil(minDef/quant)*quant:quant:ceil(maxDef/quant)*quant;
binCtrs0=binMaxs0-quant/2;
h0=hist(defs,binCtrs0);
bar(binCtrs0,h0);

[modeDef,ind]=max(h0);

for i=1:ind-1
    if h0(i)<fracMin*modeDef && i<numel(h0)
        h0(i+1)=h0(i+1)+h0(i);
        h0(i)=0;
        binMaxs0(i)=0; % delete the lower bins.
    end;
end;
for i=numel(h0):-1:ind+1
    if h0(i)<fracMin*modeDef && i>1
        h0(i-1)=h0(i-1)+h0(i);
        h0(i)=0;
        binMaxs0(i-1)=binMaxs0(i);
        binMaxs0(i)=0;
    end;
end;
q=(h0>0);
h=h0(q);
binMaxs=binMaxs0(q);
binCtrs=binMaxs-quant/2;
bar(binCtrs,h);


%%
%

nb=numel(binMaxs);
nmi=numel(allMis);
for i=1:nmi
    mi=allMis{i};
    df=mi.ctf.defocus;
    ind=find(binMaxs>=df,1);
    if numel(ind)>0
       dirName=['Info__' sprintf('%03.1f',binMaxs(ind)) '/'];
       CheckAndMakeDir(dirName,1);
       inName=allNames{i};
       [pa,nm,ex]=fileparts(inName);
       outName=[dirName nm ex];
    if doWrite
       WriteMiFile(mi,outName);
    end;
    disp(outName);
    end;
end;

return


%% Unsort
doWrite=1;
newInfoDir='Info/';
if exist(newInfoDir,'dir')
    error('Directory exists.');
end;
if doWrite
    CheckAndMakeDir(newInfoDir,1);
end;
d=dir;
for i=1:numel(d)
    if d(i).isdir && strncmp(d(i).name,'Info__',6)
        str=['cp ' d(i).name '/* ' newInfoDir];
        disp(str);
        if doWrite
            system(str);
        end;
    end;
end;