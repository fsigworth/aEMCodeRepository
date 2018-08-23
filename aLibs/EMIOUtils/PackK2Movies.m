% PackK2Movies.m
% Read 8-bit mrc movie files and convert then to packed 4-bit files.
% Add 'p4' to the end of the filename to indicate 'packed'.

% Have the user select some mrc files: boilerplate
if ~exist('fname','var') || ~exist('doBatchProcessing','var') || ~doBatchProcessing
    [fname, fpath]=uigetfile('*.mrc','Select movie files','multiselect','on');
    if isnumeric(fpath) % File selection cancelled
        return
    end;
    cd(fpath);
    if ~iscell(fname)
        fname={fname};
    end;
end;
    outputPath=fpath;
    opath=uigetdir('Where to put packed files');
    if ~isnumeric(opath) % File selection cancelled
        outputPath=AddSlash(opath);
    end;


disp(['Converting ' num2str(numel(fname)) ' files']);
for fileIndex=1:numel(fname)
    inName=fname{fileIndex};
    [pa,nm,ex]=fileparts(inName);
    outName=[outputPath nm 'p4' ex];
    disp(['Read original file ' inName]);
    [q, s]=ReadMRC(inName);
    if s.mode==32
        warning('The file is already in packed format.');
    else
        pixA=s.rez/s.nx;
        disp([' ' num2str(s.nz) ' frames.  Pixel size = ' num2str(pixA)]);
        disp(['Writing ' outName]);
        WriteMRC(q,pixA,outName,32);
    end;
    disp(' ');
end;
