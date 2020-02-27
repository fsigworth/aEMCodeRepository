% ConvertToMRCs.m
% Simple script to create mrc files from movie file formats including Tiffs.
% Scand the entire inDir for recognizable files (*.mrc;*.mrcs;*.st;*.tif;*.stk)
% Edit the following paths as necessary:

inDir='Tiffs/';
outDir='Mrcs/';

CheckAndMakeDir(outDir,1); % create output directory if needed.

d=dir(inDir);
for i=1:numel(d)
    if ~d(i).isdir % not a directory
        name=[inDir d(i).name];
        [pa,nm,ex]=fileparts(name);
            disp(['Reading ' name]);
        [m,s,ok]=ReadMovie(name);
        if ok
%             imags(BinImage(m,displayBin));  % optional image display
%             title(name,'interpreter','none');
            disp(['Writing ' outDir nm '.mrc']);
            WriteMRC(m,s.pixA,[outDir nm '.mrc']);
        end;
    end;
end;
