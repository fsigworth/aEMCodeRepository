% reMakeFlippedMaps
% For each .mrc file in the current directory, make a mirrored one with
% 'Flip' added to the filename.
d=dir;
for i=1:numel(d)
    [pa,nm,ex]=fileparts(d(i).name);
    if strcmp(ex,'.mrc')
        [m,s]=ReadMRC(d(i).name);
        outName=[nm 'Flip' ex];
        disp(outName);
        WriteMRC(MirrorX(m),s.pixA,outName);
    end;
end;