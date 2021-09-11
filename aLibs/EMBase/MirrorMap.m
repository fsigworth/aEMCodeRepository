% MirrorMap.m
% Reverse the x coordinate of the map

disp('Get the 3D map file');
[nm,pa]=uigetfile('*.mrc');
if isnumeric(pa)
    return
end;
[m0,s]=ReadMRC([pa nm]);
m1=MirrorX(m0);
[~,nm,ex]=fileparts(nm);
outName=[nm '_flip' ex];
WriteMRC(m1,s.pixA,[pa outName]);
disp([outName ' written.']);
