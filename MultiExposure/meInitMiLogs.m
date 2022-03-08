% meInitMiLogs.m
% Erase the entries in each mi.log file. This makes it easy to watch for
% incomplete processing.
infoDir='Info_C24-4/';
names=f2FindInfoFiles(infoDir);
nmi=numel(names);
nskip=0;
for i=1:nmi
    mi=ReadMiFile(names{i});
    if numel(mi.log)>0
        mi.log=cell(0);
        WriteMiFile(mi,names{i});
        disp(names{i});
    else
        nskip=nskip+1;
        if nskip>50
            nskip=0;
            fprintf(1,'\n');
        else
            fprintf(1,'.');
        end;
    end;
end;
fprintf(1,'\n');
disp('done.');
