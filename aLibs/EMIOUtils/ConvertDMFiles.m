% ConvertDMFiles.m
% Convert .dm files to .mrc and .jpg in the current directory

fileExtensions={'.dm3' '.dm4'};
d=dir;
nEntries=numel(d);
numWritten=0;
for i=3:nEntries
    name=d(i).name;
    [pa,nm,ex]=fileparts(name);  % we'll ignore the path
    if any(strcmp(ex,fileExtensions)) % We found a file of the correct type
        [m,pixA]=ReadEMFile(name);
        outName=[nm '.mrc'];
        WriteMRC(m,pixA,outName);
        disp(outName);
        n=size(m);
        md=Downsample(m,n/2);  %downsample by 2 for saving as jpeg
        jpegName=[nm,'.jpg'];
        WriteJpeg(md, jpegName);  % automatically adds the extension .jpg
        disp(jpegName);
        numWritten=numWritten+1;
    end;
end;
disp([num2str(numWritten) ' files written.']);
