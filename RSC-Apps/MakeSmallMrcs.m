% MakeSmallMrcs.m
ds=4;
d=dir;
for i=1:numel(d)
    [pa,nm,ex]=fileparts(d(i).name);
    if strcmp(ex,'.mrc')
        [m,s]=ReadMRC(d(i).name);
        nOut=size(m)/ds;
        outName=[nm 's' ex];
        disp(outName);
        WriteMRC(Downsample(m,nOut),s.pixA*ds,outName);
    end;
end;

        