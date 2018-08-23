% miFixVesicleFields.m
% Fix old version 11 mi files where vesicle.x etc. are row vectors, so that
% meMakeModelVesicles will work.

[names,pa]=uigetfile('*mi.mat','Select mi.mat files');
if isnumeric(pa)
    return
end;
cd(pa);
if ~iscell(names)
    names={names};
end;

for i=1:numel(names)
    name=names{i};
    s=load(name);
    mi=s.mi;
    fields=fieldnames(mi.vesicle);
    changed=0;
    for j=1:numel(fields)
        q=mi.vesicle.(fields{j});
        if size(q,1)==1
            changed=changed || size(q,2)>1;
            mi.vesicle.(fields{j})=q';
        end;
    end;
    if changed
        disp(['Writing ' name]);
        save(name,'mi');
    end;
end;
    