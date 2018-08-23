% LookForFiles
% Scan for newly-created files in a directory.

rexp='.*\.tif';
% rexp='.*\.mrc'
names={};
a=struct;
a.FileSize=0;
newName='';
ok=1;
while ok
    d=dir;
    found=0;
    for i=1:numel(d)
        q=regexp(d(i).name,rexp);
        if ~d(i).isdir && numel(q)>0
            if ~any(strcmp(d(i).name,names))  % a new name
                found=1;
                names{end+1}=d(i).name;
                newName=d(i).name
                j=i;
            end;
        end;
    end;
    if numel(newName)>0
        try
%         s=warning('off','all');
        a=imfinfo(newName);
%         warning(s);
        disp([a(1).FileSize d(j).bytes]);
        end;
%         f=fopen(newName)
%         fclose(f);
    end;
        pause(0.1);
end;
