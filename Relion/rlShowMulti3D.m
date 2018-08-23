% rlShowMulti3D.m

[names,path]=uigetfile('*.mrc*','multiselect','on');
if isnumeric(path)
    return
end;
cd(path);

if ~isa(names,'cell')
    names={names};
end;
%%
maps=[];
for i=1:numel(names)
    disp(names{i});
    maps(:,:,:,i)=ReadMRC(names{i});
end;

ShowSections(maps,[],45); % not using opts because ShowSections can't make labels.
