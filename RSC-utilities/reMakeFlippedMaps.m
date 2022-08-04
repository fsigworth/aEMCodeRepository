function reMakeFlippedMaps(names)
% function reMakeFlippedMaps(names)
% flips 3D maps along X. 
% - If no argument given, it creates a new file for
% every .mrc file in the current directory.
% - If names is a character string, flips just that name.

% For each .mrc file in the current directory, make a mirrored one with
% 'Flip' added to the filename.
if nargin<1 | numel(names)<1
   d=dir;
    for i=1:numel(d)
        names{i}=d(i).name;
    end;
elseif isa(names,"char")
    names={names};
end;
nd=numel(names);

for i=1:nd
    [pa,nm,ex]=fileparts(names{i});
    if strcmp(ex,'.mrc')
        [m,s]=ReadMRC(names{i});
        outName=[AddSlash(pa) nm 'Flip' ex];
        disp(outName);
        WriteMRC(MirrorX(m),s.pixA,outName);
    end;
end;
