% rlCTFsToJpeg.m

d=dir;
outDir='ctf-jpg/';

for i=1:numel(d)
    name=d(i).name;
    if ~d(i).isdir && numel(regexp(name,'\w+\.ctf'))>0
        % has .ctf in the name
        disp(name);
        CheckAndMakeDir(outDir,1);
        m=ReadMRC(name);
        [pa,nm,ex]=fileparts(name);
        outName=[outDir nm '.jpg'];
        WriteJpeg(m',outName,0);
    end;
end;
