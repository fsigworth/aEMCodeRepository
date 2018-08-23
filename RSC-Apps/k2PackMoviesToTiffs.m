% k2PackMoviesToTiffs.m
% Read 8-bit mrc movie files and convert then to lzw-compressed tiff files.
% k2 movies get compressed to about 1/3 original size.
doDeleteOriginal=1
doUpdateMi=1

if ~exist('gNames','var') || ~exist('gBatchProcessing','var') || ~gBatchProcessing
    [gNames, pa]=uigetfile('*mi.mat','Select mi files','multiselect','on');
    if isnumeric(pa) % File selection cancelled
        return
    end;
    [rootPath, gInfoPath]=ParsePath(pa);
    cd(rootPath);
    if ischar(gNames)
        gNames={gNames};
    end;
end;

tic
disp(['Converting ' num2str(numel(gNames)) ' files']);
for fileIndex=1:numel(gNames)
    miName=gNames{fileIndex};
    load([gInfoPath miName]);
    disp([num2str(fileIndex) ':  ' miName]);
    dataPath=mi.moviePath;
    inName=mi.movieFilename;
    [pa,nm,ex]=fileparts(inName);
    if strcmp(ex,'.mrc')  % needs conversion
        if strcmp(nm(end-1:end),'p4')  % remove the p4 suffix
            nm(end-1:end)=[];
        end;
        outName=[nm '.tif'];
        disp(['Read original file ' dataPath inName]);
        [q, s]=ReadMRC([dataPath inName]);
        pixA=s.pixA;
        disp([' ' num2str(s.nz) ' frames.  Pixel size = ' num2str(pixA)]);
        disp(['Writing ' dataPath outName]);
        WriteTiffStack(q,pixA,[dataPath outName]);
        if doDeleteOriginal
            delete([dataPath inName]);
        end;
        if doUpdateMi
            mi.movieFilename=outName;
            mi.basePath=rootPath;
            save([gInfoPath miName],'mi');
        end;
        disp(' ');
    end;
end;
toc