% SearchRedundantMis.m

d=dir('Info/');
names={};
strs={};
j=0;
for i=1:numel(d)
    name=d(i).name;
    if strndcmp(name,'mi.txt')
        tstring=name(11:23); % part of the name we'll use for match
        j=j+1;
        names{j}=name;
        strs{j}=tstring;
    end;
end;
[sstrs,inds]=sort(strs);
for i=2:j
    if strcmp(sstrs{i},sstrs{i-1})
        name=names{inds(i)};
        infoCmd=['!mv Info/' name ' InfoExtra/'];
%         eval(infoCmd);
        baseName=name(1:end-6);
        mergeCmd=['!mv Merged/' baseName '* MergeExtra/'];
%         eval(mergeCmd);
        disp([infoCmd '  ' mergeCmd]);
    end;
end;
% for i=1:j
%     disp(names{inds(i)});
% end;


%%
d=dir('InfoExtra/');
baseNames={};
j=0;
for i=1:numel(d)
    name=d(i).name;
    if strndcmp(name,'mi.txt')
        j=j+1;
        baseNames{j}=name(1:end-6);
    end;
end;
for i=1:j
%     disp(baseNames{i});
            mgCmd=['!mv Micrograph/' baseNames{i} '* MicrographExtra/'];
            eval(mgCmd);
disp(mgCmd);
end;