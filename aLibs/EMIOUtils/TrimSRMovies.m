% TrimSRMovies.m
% Trim super-res movies to [7676 7420] size.
% Reads every *.tif file in the current directory.  If it has the correct
% original size, crops each frame and writes out the stack to an .mrcs file
% in 8-bit format.

origSize=[7680 7424];
trimSize=[7676 7420];

d=dir;
for i=1:numel(d)
    name=d(i).name;
    [pa,nm,ex]=fileparts(name);
    if strcmp(ex,'.tif') % it's a tif, presumably movie file.
        disp(['Reading ' name]);
        m=ReadMovie(name);
        sz=size(m);
        disp([' size: ' num2str(sz)]);
        if all(sz(1:2)==origSize)   % Check that it's the original size
            mcr=Crop(m,trimSize,1); % trim each frame in the stack
            outName=[nm '.mrcs'];
            disp([' Writing ' outName]);
            WriteMRC(mcr,0,outName,0); % pixA=0; mode 0 means uint8
        else
            disp(' -wrong size, skipped');
        end;
    end;
end;
