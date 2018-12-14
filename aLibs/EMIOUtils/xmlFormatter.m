% xmlFormatter.m
% Split an xml file into lines.

[nm,pa]=uigetfile('*.xml');
if isnumeric(pa)
    return
end;
cd(pa);
f1=fopen(nm);
a=fread(f1,inf)';
fclose(f1);

a=char(a);
ptrs1=strfind(a,'<');
ptrs2=strfind(a,'</');
ptrs3=strfind(a,'>');
p21=1;
for i=1:numel(ptrs2)
    p2=ptrs2(i); % current end begins with </...
    q22=find(ptrs1>p2,1); % next opening
    if numel(q22)>0
        p22=ptrs1(q22); % next opening
    else
        p22=0;
    end;
    if any(ptrs2==p22) % that opening is actually an end
        continue;
    end;
    q1=find(ptrs1>p21,1); % first opening following last end;
    q3=find(ptrs3>p2,1);
    p1=ptrs1(q1);
    p3=ptrs3(q3);
    disp(a(p1:p3));
    p21=p2;
end;

