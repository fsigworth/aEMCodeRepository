% rlRescaleMC2Micrographs.m
% Pad to our favorite size and scale (if needed) DW micrographs created by
% MotionCor2.
scale=4;

d=dir;
for i=1:numel(d)
    name=d(i).name;
    [pa,nm,ex]=fileparts(name);
    if strcmp(ex,'.mrc') && nm(end)~='c'
        [m,s]=ReadMRC(name);
        m=m*scale;
        me=mean(double(m(:)));
        sz=size(m);
        newSz=NextNiceNumber(sz,5,4);
        m1=Crop(m,newSz,0,me);
        binFactor=max(1,round(newSz(1)/500));
        imags(BinImage(m1,binFactor));
        title(name,'interpreter','none');
        drawnow;
        outName=[nm 'c' ex];
        WriteMRC(m1,s.pixA,outName);
        disp([outName ' written. ' num2str(me,6)]);
        return
    end;
end;
